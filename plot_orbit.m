%==========================================================================
% Interplanetary Trajectory Optimization for Planetary Defense Mission
% (Asteroid Kinetic Delfection)
% 
% Lee Kin Thong 
% Sept 22 2025
%==========================================================================
% You are free to use and modify the code, but you MUST cite the following
% papers:
%
% Lee, Kinthong, Zhengqing Fang, and Zhaokui Wang. "Investigation of the 
% incremental benefits of eccentric collisions in kinetic deflection of 
% potentially hazardous asteroids." Icarus 425 (2025): 116312.
%
% Feels free to contact me for any inquiry or cooperation!
% ktlee3819@gmail.com
%==========================================================================
% This code will plot 3 figures and 1 animation as used in the paper:
% Figures:
% 1. Low-thrust acceleration in-plane & out-of-plane angle
% 2. Transfer trajectory, both detailed(large, with low-thrust) and small
%    (small without low-thrust)
% 3. Intercept moment to show interception angle
%
% Animation:
% 1. Intercept animation
%--------------------------------------------------------------------------



clear

%--------------------------------------------------------------------------
% ---------------------------  Add Path  ----------------------------------
%--------------------------------------------------------------------------
currentDir = fileparts(mfilename('fullpath'));
addpath(genpath(currentDir))

% different path separators for Window(\), MacOS(/), and Linux(/) 
if ispc
    pathTosaveLowthrust = [fullfile(currentDir, 'output_result','figure','lowthrust'),'\'];
    pathTosaveTransferTrajectory = [fullfile(currentDir, 'output_result','figure','transfer_trajectory'),'\'];
    pathTosaveDetailedTrajectory = [fullfile(currentDir, 'output_result','figure','detailed_trajectory'),'\'];
    pathTosaveInterception = [fullfile(currentDir, 'output_result','figure','interception'),'\'];
    pathTosaveAnimation = [fullfile(currentDir, 'output_result','animation','transfer_trajectory'),'\'];
else
    pathTosaveLowthrust = [fullfile(currentDir, 'output_result','figure','lowthrust'),'/'];
    pathTosaveTransferTrajectory = [fullfile(currentDir, 'output_result','figure','transfer_trajectory'),'/'];
    pathTosaveDetailedTrajectory = [fullfile(currentDir, 'output_result','figure','detailed_trajectory'),'/'];
    pathTosaveInterception = [fullfile(currentDir, 'output_result','figure','interception'),'/'];
    pathTosaveAnimation = [fullfile(currentDir, 'output_result','animation','transfer_trajectory'),'/'];
end

% add path to kernel  
% Specify the directory where your .bsp files are located
pathTokernel = fullfile(currentDir, 'code', 'kernel');
%--------------------------------------------------------------------------





% -------------------------------------------------------------------------
%-----------------------------Load Kernels---------------------------------
% -------------------------------------------------------------------------
% Load Kernels not in parfor loop
% List all .bsp files in the directory
bspFiles = dir(fullfile(pathTokernel, '*.bsp'));

% List all .txt files in the directory
txtFiles = dir(fullfile(pathTokernel, '*.txt'));

% Combine the .bsp and .txt files into a single array
allFiles = [bspFiles; txtFiles];

% Iterate over all files to load them
for k = 1:length(allFiles)
    % Construct the full path for each file
    filePath = fullfile(pathTokernel, allFiles(k).name);

    % Load the file using cspice_furnsh
    cspice_furnsh(filePath);
end
% -------------------------------------------------------------------------




%--------------------------------------------------------------------------
%------------------- Read data from PHA_table.xlsx ------------------------
%--------------------------------------------------------------------------
% Read PHA_table.xlsx
filename = 'PHA_table.xlsx';
data = readtable(filename, 'PreserveVariableNames', true);

% Extract PHA name and ID code
PHA_name = string(data.Object);
PHA_bsp = string(data.BSP_file_name); % example: 2099942: Apophis

% Extract year, month, and day from the 'Close_Approach_CA_Date' column
dateStrings = data.Close_Approach_CA_Date;
datePattern = '(\d{4})-(\w{3})-(\d{2})';
tokens = regexp(dateStrings, datePattern, 'tokens');

% Convert tokens to a matrix (each row corresponds to a match)
tokensMatrix = vertcat(tokens{:});
tokensMatrix = vertcat(tokensMatrix{:});

% Extract and convert each component
year_CA = str2double(tokensMatrix(:,1)); 
date_CA = str2double(tokensMatrix(:,3));

% Map the month abbreviations to month numbers and full names
monthAbbreviations = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
monthNumbers = 1:12;
monthFullNames = {'January', 'February', 'March', 'April', 'May', 'June', ...
                  'July', 'August', 'September', 'October', 'November', 'December'};

% Create containers.Map to map month abbreviations to month numbers and full names
monthNumMap = containers.Map(monthAbbreviations, monthNumbers);
monthNameMap = containers.Map(monthAbbreviations, monthFullNames);

% Convert month abbreviations to numbers and full names
month_CA_num = cell2mat(values(monthNumMap, tokensMatrix(:,2)));
month_CA_string = values(monthNameMap, tokensMatrix(:,2));

% Extract position result
% Extract the first 32 rows of 'BestPosition_COG' as a cell array of strings
A_COG = data.BestPosition_COG(1:32);
A_BIP = data.BestPosition_BIP(1:32);

% Initialize the matrix
A_COG1 = zeros(32, 17);
A_BIP1 = zeros(32, 17);

% Process each row and convert to numeric format
for i = 1:32
    str_COG = A_COG{i}; % Get the string
    str_BIP = A_BIP{i}; % Get the string
    
    % Remove any square brackets if present
    str_COG = erase(str_COG, {'[', ']'});
    str_BIP = erase(str_BIP, {'[', ']'});
    
    % Convert the string of numbers into a numeric row vector
    % the 18th column in NaN
    A_COG1(i, :) = str2double(split(str_COG, ','))';
    A_BIP1(i, :) = str2double(split(str_BIP, ','))';
end
Position_COG = A_COG1(1:32,1:17);
Position_BIP = A_BIP1(1:32,1:17);
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% --------------------------- Parameters ----------------------------------
% -------------------------------------------------------------------------
% Modify as your need

% Units
AU = 149597870691; % m/AU
TU = 5022548.032; % s/TU

% Parameters of collision model
problem.mass_PHA = 2.7e10; % kg
problem.beta = 3.61;

% Parameters for rocket
problem.c_impulse = 3.09; % Exhaust velocity km/s, I = 315s
problem.c_lowthrust = 29.78  ; % Exhaust velocity km/s
problem.m0 = 10000; % kg
problem.n0 = 18e-6 * 9.81 ; % thrust-to-mass ratio at t0, m/s^2
problem.I_rocket = 3035; % s

% PHA 
target = data.BSP_file_name;
% SOI date
date_SOI = data.Date_SOI;


% -------------------------------------------------------------------------
% ---------------------------- Propagation --------------------------------
% -------------------------------------------------------------------------
% PHA number, as label in PHA_table.xlsx
for PHA = 1:1
    for is_BIP = 0 : 1
        % Problem defination
        problem.is_BIP = is_BIP;
        problem.Target = char(PHA_bsp(PHA)); % Target, 2099942 represents Apophis
        problem.date_SOI = data.Date_SOI(PHA); % Mjday
        problem.date_CA = Mjday(year_CA(PHA),month_CA_num(PHA),date_CA(PHA));
        problem.CA_distance = data.CADistanceMinimum_km_(PHA) * 1000; % m
        
        % Parameters of Particle
        if is_BIP
            Particle_Params.thetaIN_coeff = Position_BIP(PHA,1:6);
            Particle_Params.thetaOUT_coeff = Position_BIP(PHA,7:12);
            Particle_Params.v_inf_abs = Position_BIP(PHA,13);
            Particle_Params.RA = Position_BIP(PHA,14);
            Particle_Params.DEC = Position_BIP(PHA,15);
            Particle_Params.launch_date = Position_BIP(PHA,16);
            Particle_Params.impact_date = Position_BIP(PHA,17);
        else
            Particle_Params.thetaIN_coeff = Position_COG(PHA,1:6);
            Particle_Params.thetaOUT_coeff = Position_COG(PHA,7:12);
            Particle_Params.v_inf_abs = Position_COG(PHA,13);
            Particle_Params.RA = Position_COG(PHA,14);
            Particle_Params.DEC = Position_COG(PHA,15);
            Particle_Params.launch_date = Position_COG(PHA,16);
            Particle_Params.impact_date = Position_COG(PHA,17);
        end
% -------------------------------------------------------------------------


        
        
% ------------------ Propagation Low Thrust to Impact ---------------------
        % Obtain Initial State
        % Obtain initial conditions from SPICE
        % Convert date from MJday to ET (Ephemeris Time)

        % Launch Date
        [years_launch,months_launch, days_launch, fds_launch] = iauJd2cal( 2400000.5, Particle_Params.launch_date);
        [hours_launch, minutes_launch, secs_launch] = fd_to_hms(fds_launch);
        [years_launch, months_launch, days_launch, hours_launch, minutes_launch, secs_launch] = fix_seconds(years_launch, months_launch, days_launch, hours_launch, minutes_launch, secs_launch);
        secs_launch = round(secs_launch,3);
        months_launch_string = monthToString(months_launch);
        str = sprintf(' %s %g , %g %g:%g:%g', months_launch_string, days_launch, years_launch, hours_launch, minutes_launch, secs_launch);
        et_launch = cspice_str2et(str);
        [Y_launch, ~] = cspice_spkezr('EARTH', et_launch, 'J2000', 'NONE', 'SUN');

        % Calculate the impulse delta_v at the launch date
        % Initial State of spacecarft
        r = [Y_launch(1);Y_launch(2);Y_launch(3)];
        v = [Y_launch(4);Y_launch(5);Y_launch(6)];
        
        % Construct the rotation matrix to transform to heliocentric frame
        % Unit vectors of the local orbital frame
        r_hat = r / norm(r); % Radial unit vector
        h = cross(r, v);     % Specific angular momentum vector
        h_hat = h / norm(h); % Normal unit vector (orbital plane normal)
        t_hat = cross(h_hat, r_hat); % Tangential unit vector
        
        % Rotation matrix (orbital frame to heliocentric frame)
        R = [r_hat, t_hat, h_hat]; % Each column is a unit vector
        
        % Transform acceleration to heliocentric frame
        delta_v = R * [cos(Particle_Params.DEC)*cos(Particle_Params.RA);cos(Particle_Params.DEC)*sin(Particle_Params.RA);sin(Particle_Params.DEC)]; % Transform orbital to heliocentric frame
        
        % Add v_inf to initial velocity
        V_impulse = delta_v .* Particle_Params.v_inf_abs;
        Y_launch(4:6) = Y_launch(4:6) + V_impulse;
        Y_launch = Y_launch.*1000; % m, m/s
        Y_launch = Y_launch ./ AU; % AU, AU/s
        Y_launch(4:6) = Y_launch(4:6) .* TU; % AU/TU
        

        % Impact Date
        [years_impact,months_impact, days_impact, fds_impact] = iauJd2cal( 2400000.5, Particle_Params.impact_date);
        [hours_impact, minutes_impact, secs_impact] = fd_to_hms(fds_impact);
        [years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact] = fix_seconds(years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact);
        secs_impact = round(secs_impact,3);
        months_impact_string = monthToString(months_impact);
        str = sprintf(' %s %g , %g %g:%g:%g', months_impact_string, days_impact, years_impact, hours_impact, minutes_impact, secs_impact);
        et_impact = cspice_str2et(str);
        [Y_PHA_impact, ~] = cspice_spkezr(problem.Target, et_impact, 'J2000', 'NONE', 'SUN');
        Y_PHA_impact = Y_PHA_impact.*1000; % m, m/s
        Y_PHA_impact = Y_PHA_impact ./ AU; % AU, AU/s
        Y_PHA_impact(4:6) = Y_PHA_impact(4:6) .* TU; % AU/TU
        
    
        
        % Make Sure Impact date is larger than Launch Date
        % Calculate low thrust transfer tracjectory to obtain v_sc at impact &
        % interception distance
        tspan = [0, (et_impact - et_launch)/TU];

        % Starts Propagation
        Y_spacecraft_launch_impact = Ephemeris_lowthrust_heliocentric(Y_launch, tspan,Particle_Params,problem);
    
        % Spacecraft State at impact moment
        constrain_for_collision = norm(Y_spacecraft_launch_impact(end,2:4)' - Y_PHA_impact(1:3));
        Interception = constrain_for_collision * AU; % m
        InterceptionAngle = acosd( dot( Y_spacecraft_launch_impact(end,5:7)',Y_PHA_impact(4:6) ) / norm(Y_spacecraft_launch_impact(end,5:7)') / norm(Y_PHA_impact(4:6)) ); % deg
    
    
% ------------------------------- Impact ----------------------------------
        % Calculate remaining mass
        mass_after_depart_Earth = problem.m0 * exp( - Particle_Params.v_inf_abs / problem.c_impulse ); % kg
        mass_at_impact = mass_after_depart_Earth * exp( - problem.n0 / (problem.c_lowthrust*1000) * (et_impact - et_launch) ); % kg
    
        % Calculate delta_v after impact
        % Read 3D models
        [vertices, faces] = readObj('Apophis_Model.obj');
        % Collision Model
        if problem.is_BIP 
            % BIP stategy
            [delta_v,~,~] = get_delta_v_from_momentum(Y_spacecraft_launch_impact(end,5:7)', Y_PHA_impact(4:6), mass_at_impact ,problem,vertices,faces);
            delta_v_to_PHA_at_impact = delta_v' .* AU ./ TU; % m, m/s, switch from 1x3 to 3x1
            save_type = 'BIP';
            data.Launch_BIP(PHA) = datetime(years_launch,months_launch, days_launch,hours_launch, minutes_launch, secs_launch);
            data.Impact_BIP(PHA) = datetime(years_impact,months_impact, days_impact,hours_impact, minutes_impact, secs_impact);
            data.Flight_Duration_BIP(PHA) = (et_impact - et_launch)/86400;
            data.Impulse_BIP(PHA) = Particle_Params.v_inf_abs;
            data.Interception_Angle_BIP(PHA) = InterceptionAngle;
            data.Interception_Error_BIP(PHA) = Interception;
            data.Impact_Mass_BIP(PHA) = mass_at_impact;
        else
            % COG strategy
            [~,delta_v,~] = get_delta_v_from_momentum(Y_spacecraft_launch_impact(end,5:7)', Y_PHA_impact(4:6), mass_at_impact ,problem,vertices,faces);
            delta_v_to_PHA_at_impact = delta_v .* AU ./ TU; % m, m/s, NO NEED switch, it is already 3x1
            save_type = 'COG';
            data.Launch_COG(PHA) = datetime(years_launch,months_launch, days_launch,hours_launch, minutes_launch, secs_launch);
            data.Impact_COG(PHA) = datetime(years_impact,months_impact, days_impact,hours_impact, minutes_impact, secs_impact);
            data.Flight_Duration_COG(PHA) = (et_impact - et_launch)/86400;
            data.Impulse_COG(PHA) = Particle_Params.v_inf_abs;
            data.Interception_Angle_COG(PHA) = InterceptionAngle;
            data.Interception_Error_COG(PHA) = Interception;
            data.Impact_Mass_COG(PHA) = mass_at_impact;
        end
% -------------------------------------------------------------------------
    
    



% -------------- Propagation from impact date to CA date ------------------
        % To calculate the deflection distance,
        % The proppagation is divided into two step:
        % 1. from Impact date to Earth SOI : in heliocentric, a step size of 1 day (86400secs)
        % 2. from Earth SOI to Ca date + 1 day : in geocentric, a step size of 1 seconds

        %  Propagation from impact date to closest-approach date
        DeflectionDistance = Propagation_impact_CA(delta_v_to_PHA_at_impact,Particle_Params,problem); % m

        % Write to data
        if problem.is_BIP 
            data.BestDeflectionDistance_BIP(PHA) = DeflectionDistance;
        else
            data.BestDeflectionDistance_COG(PHA) = DeflectionDistance;
        end


        % Write to table
        writetable(data,'temporary_result1.xlsx')

            

% -------------------------- Store date for plot --------------------------
        % Obtain position of Earth and PHA
        % launch to impact time
        et = et_launch +  Y_spacecraft_launch_impact(:,1)' .* TU;
        % Earth state, launch to impact
        [Y_Earth_launch_impact, ~] = cspice_spkezr('Earth', et, 'J2000', 'NONE', 'SUN');
        Y_Earth_launch_impact = Y_Earth_launch_impact .* 1000 ./ AU;
        % PHA state, launch to impact
        [Y_PHA_launch_impact, ~] = cspice_spkezr(problem.Target, et, 'J2000', 'NONE', 'SUN');
        Y_PHA_launch_impact = Y_PHA_launch_impact .* 1000 ./ AU;

        % Orbit period
        [a,~,~,~,~,~] = input_r_v(Y_PHA_launch_impact(1:3,1) .* AU ,Y_PHA_launch_impact(4:6,1) .* AU, 132712440041279422464);
        T = 2 * pi * sqrt( a^3 / 132712440041279422464 );
        % PHA's heliocentric orbit
        et = et_launch : 86400 : et_launch + T + 5*86400 ;
        [Y_PHA_orbit_period, ~] = cspice_spkezr(problem.Target, et, 'J2000', 'NONE', 'SUN');
        Y_PHA_orbit_period = Y_PHA_orbit_period .* 1000 ./ AU;
        % Earth's heliocentric orbit
        et = et_launch : 86400 : et_launch + 86400*366 ;
        [Y_Earth_orbit_period, ~] = cspice_spkezr('Earth', et, 'J2000', 'NONE', 'SUN');
        Y_Earth_orbit_period = Y_Earth_orbit_period .* 1000 ./ AU;
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
% --------------- Plot low thrust accerelation angle ----------------------
% -------------------------------------------------------------------------
        % Thrust components in spacecraft orbit coordinate
        syms t
        % In-plane thrust pointing angle (thetaIN)
        thetaIN(t) = dot(Particle_Params.thetaIN_coeff , [1 ,t ,t^2 ,t^3, t^4, t^5] );
        % Out-of-plane thrust pointing angle (thetaOUT)
        thetaOUT(t) = dot(Particle_Params.thetaOUT_coeff , [1 ,t ,t^2 ,t^3, t^4, t^5] );

         a_low_thrust = zeros(3,length(Y_spacecraft_launch_impact));
         for i = 1 : length(Y_spacecraft_launch_impact)
             time = Y_spacecraft_launch_impact(i,1);
            a_r = cos(thetaOUT(time)) * cos(thetaIN(time)); % Radial acceleration
            a_t = cos(thetaOUT(time)) * sin(thetaIN(time)); % Tangential acceleration
            a_n = sin(thetaOUT(time));             % Normal acceleration

            % Acceleration vector in the local orbital frame
            a_orbital = [a_r; a_t; a_n];

            % Construct the rotation matrix to transform to heliocentric frame
            % Unit vectors of the local orbital frame
            r = [Y_spacecraft_launch_impact(i,2);Y_spacecraft_launch_impact(i,3);Y_spacecraft_launch_impact(i,4)];
            v = [Y_spacecraft_launch_impact(i,5);Y_spacecraft_launch_impact(i,6);Y_spacecraft_launch_impact(i,7)];
            r_hat = r / norm(r); % Radial unit vector
            h = cross(r, v);     % Specific angular momentum vector
            h_hat = h / norm(h); % Normal unit vector (orbital plane normal)
            t_hat = cross(h_hat, r_hat); % Tangential unit vector

            % Rotation matrix (orbital frame to heliocentric frame)
            R = [r_hat, t_hat, h_hat]; % Each column is a unit vector

            % Transform acceleration to heliocentric frame
            a_low_thrust(:,i) = R * a_orbital; % Transform orbital to heliocentric frame
         end

        % Plot
        fig1 = figure(1);

        % In-plane angle
        subplot(2,1,1)
        % Convert to deg
        thetaIN(t) = mod(thetaIN(t) * 180/pi + 180 ,360) - 180;
        fplot(thetaIN,[0,Y_spacecraft_launch_impact(end,1)])
        title('In-plane thrust pointing angles \beta')
        xlabel('Duration (TU)')
        ylabel('Angle (deg)')

        % Out-of-plane angle
        subplot(2,1,2)
        thetaOUT(t) = mod(thetaOUT(t) * 180/pi + 180 ,360) -180;
        fplot(thetaOUT,[0,Y_spacecraft_launch_impact(end,1)])
        title('Out-of-plane thrust pointing angles \gamma')
        xlabel('Duration (TU)')
        ylabel('Angle (deg)')


        fig1.Units = 'inches';
        fig1.Position = [0, 1, 7, 7]; % [left, bottom, width, height]

        % Export the png file 
        exportgraphics(fig1, [pathTosaveLowthrust,num2str(PHA),save_type,'_low_thrust_angle.png']);
        close(fig1)
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% -------------------- Plot Transfer Trajectory  --------------------------
% -------------------------------------------------------------------------
        % Figure 2: Small figure, without low thrust acceleration
        fig2 = figure(2);
        % Figure 3: Large figure, with low thrust acceleration
        fig3 = figure(3);

        % For labeling, fXXYY:
        % XX: 02: Small figure (fig2), 03: Large figure (fig3)
        % YY: number of elements plotted


        % Impulse acceleration unit vector
        V_impulse =  V_impulse ./ norm( V_impulse);


        % The code will go for COG first (else loop), then BIP (if loop)
        if is_BIP
            % Spacecraft
            figure(2)
            f0207 = plot3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4),'LineWidth',2,'Color', [1, 0.5, 0]);
            hold on

            figure(3)
            f0309 = plot3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4),'LineWidth',2,'Color', [1, 0.5, 0]);
            hold on

            % Larger detailed figure
            figure(3)
            % Impulse at begining
            f0310 = quiver3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4), V_impulse(1), V_impulse(2), V_impulse(3), 'b', 'AutoScaleFactor', 0.2);
            % Accerelation
            f0311 = quiver3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4), a_low_thrust(1,:)', a_low_thrust(2,:)', a_low_thrust(3,:)','color', [1, 0.5, 0], 'AutoScaleFactor', 0.5);
            hold on
            % Launch Point
            f0312 = scatter3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4),200, [1, 0.5, 0],'filled');
            hold on
            % Impact Point
            f0313 = scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),500, [1, 0.5, 0],'pentagram','filled');
            hold on

            title(sprintf('%d. %s', PHA, PHA_name(PHA)));
            xlabel('X (AU)')
            ylabel('Y (AU)')
            zlabel('Z (AU)') 

            view(-45,45)
            xL = xlim;  
            yL = ylim;  
            xlim([min([xL(1), yL(1)]), max([xL(2), yL(2)])]);
            ylim([min([xL(1), yL(1)]), max([xL(2), yL(2)])]);
            fig3.Units = 'inches';

            fig3.Position = [1, 0, 16, 11]; % [left, bottom, width, height]
            pause(0.5)
            lg = legend([f0301, f0302,f0303,f0306,f0309,f0311,f0305], {'Earth Orbit','PHA Orbit', 'COG trajectory','COG Low-Thrust Accerelation', 'BIP trajectory', 'BIP Low-Thrust Accerelation', 'Impulse'});
            lg.Units = 'inches';
            lg.Position = [13,6,0,0];
            



            % -------------------Small Picture--------------------
            % To better fit the paper, the images have been restricted to a small size, but they can be adjusted as needed.
            % To plot small picture to suit into paper
            figure(2)
            % Launch Point
            f0208 = scatter3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4),50, [1, 0.5, 0],'filled');
            hold on
            % Impact Point
            f0209 = scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),50, [1, 0.5, 0],'pentagram','filled');
            hold on

            title(sprintf('%d. %s', PHA, PHA_name(PHA)));
            xlabel('X (AU)')
            ylabel('Y (AU)')
            zlabel('Z (AU)') 

            view(-45,45)
            xL = xlim;  
            yL = ylim;  
            xlim([min([xL(1), yL(1)]), max([xL(2), yL(2)])]);
            ylim([min([xL(1), yL(1)]), max([xL(2), yL(2)])]);

            fig2.Units = 'inches';
            fig2.Position = [1, 1, 3, 3]; % [left, bottom, width, height]

            % Give time for matlab to render it, if not the size might
            % mismatch
            pause(0.5)
            

            % Export graphics
            exportgraphics(fig2, [pathTosaveTransferTrajectory,num2str(PHA),'_Transfer_Orbit.png']);
            close
            exportgraphics(fig3, [pathTosaveDetailedTrajectory,num2str(PHA),'_Transfer_Orbit.png']);
            close



        else

            % ------------------- Small figure ----------------------------
            figure(2)
            % Earth's Orbit
            f0201 = plot3(Y_Earth_orbit_period(1,:),Y_Earth_orbit_period(2,:),Y_Earth_orbit_period(3,:),'LineWidth',2,'Color','Blue');
            hold on
            % PHA's Orbit
            f0202 = plot3(Y_PHA_orbit_period(1,:),Y_PHA_orbit_period(2,:),Y_PHA_orbit_period(3,:),'LineWidth',2,'Color','Green');
            hold on
            % Spacecraft trajectory
            f0203 = plot3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4),'LineWidth',2,'Color','Red');
            hold on
            % Sun position
            f0204 = scatter3(0,0,0,'Red');
            hold on
            % Launch Point
            f0205 = scatter3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4),50,'red','filled');
            hold on
            % Impact Point
            f0506 = scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),50,'pentagram','red','filled');
            hold on

            % ------------------- Large figure ----------------------------
            figure(3)
            % Earth's Orbit
            f0301 = plot3(Y_Earth_orbit_period(1,:),Y_Earth_orbit_period(2,:),Y_Earth_orbit_period(3,:),'LineWidth',2,'Color','Blue');
            hold on
            % PHA's Orbit
            f0302 = plot3(Y_PHA_orbit_period(1,:),Y_PHA_orbit_period(2,:),Y_PHA_orbit_period(3,:),'LineWidth',2,'Color','Green');
            hold on
            % Spacecraft trajectory
            f0303 = plot3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4),'LineWidth',2,'Color','Red');
            hold on
            % Sun position
            f0304 = scatter3(0,0,0,'Red');
            hold on

            % Impulse at begining
            f0305 = quiver3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4), V_impulse(1), V_impulse(2), V_impulse(3), 'b', 'AutoScaleFactor', 0.2);
            hold on
            % Low-thrust Accerelation
            f0306 = quiver3(Y_spacecraft_launch_impact(:,2),Y_spacecraft_launch_impact(:,3),Y_spacecraft_launch_impact(:,4), a_low_thrust(1,:)', a_low_thrust(2,:)', a_low_thrust(3,:)', 'r', 'AutoScaleFactor', 0.5);
            hold on
            % Launch Point
            f0307 = scatter3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4),200,'red','filled');
            hold on
            % Impact Point
            f0308 = scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),500,'pentagram','red','filled');
            hold on


        end
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% ------------------------- Plot Interception -----------------------------
% -------------------------------------------------------------------------
        fig4 = figure(4);
        % PHA
        plot3(Y_PHA_launch_impact(1,end-1:end),Y_PHA_launch_impact(2,end-1:end),Y_PHA_launch_impact(3,end-1:end),'LineWidth',2,'Color','Green')
        hold on
        if is_BIP
            % Spacecraft
            plot3(Y_spacecraft_launch_impact(end-1:end,2),Y_spacecraft_launch_impact(end-1:end,3),Y_spacecraft_launch_impact(end-1:end,4),'LineWidth',2,'Color',[1, 0.5, 0])
            hold on
            % Impact
            scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),500,[1 0.5 0],'pentagram','filled')
            hold on
        else
            % Spacecraft
            plot3(Y_spacecraft_launch_impact(end-1:end,2),Y_spacecraft_launch_impact(end-1:end,3),Y_spacecraft_launch_impact(end-1:end,4),'LineWidth',2,'Color','Red')
            hold on
            % Impact
            scatter3(Y_spacecraft_launch_impact(end,2),Y_spacecraft_launch_impact(end,3),Y_spacecraft_launch_impact(end,4),500,'pentagram','red','filled')
            hold on
        end



        title('Interception')
        xlabel('X (AU)')
        ylabel('Y (AU)')
        zlabel('Z (AU)')


        % View from cross product of both vector
        v1 = [Y_spacecraft_launch_impact(end,5);Y_spacecraft_launch_impact(end,6);Y_spacecraft_launch_impact(end,7)];
        v2 = [Y_PHA_launch_impact(4,end);Y_PHA_launch_impact(5,end);Y_PHA_launch_impact(6,end)];
        normalVector = cross(v1, v2);
        view(normalVector)
        % To better fit the paper, the images have been restricted to a small size, but they can be adjusted as needed.
        fig4.Units = 'inches';
        % lg = legend('PHA Orbit','Spacecraft Trajectory');
        fig4.Position = [3, 3, 6, 5]; % [left, bottom, width, height]
        % lg.Units = 'inches';
        % lg.Position = [3,4,0,0];

        % Export the png file to specific file
        exportgraphics(fig4, [pathTosaveInterception,num2str(PHA),save_type,'_Interception.png']);
        close
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% ---------------------- Plot transfer animation --------------------------
% -------------------------------------------------------------------------
        fig = figure(5);
        % To better fit the paper, the images have been restricted to a small size, but they can be adjusted as needed.
        fig.Units = 'inches';
        fig.Position = [0, 1, 15, 11]; % [left, bottom, width, height]

        % plot first point for legend
        % Orbit
        plot3(Y_Earth_launch_impact(1,1),Y_Earth_launch_impact(2,1),Y_Earth_launch_impact(3,1),'LineWidth',2,'Color','Blue');
        hold on
        plot3(Y_PHA_launch_impact(1,1),Y_PHA_launch_impact(2,1),Y_PHA_launch_impact(3,1),'LineWidth',2,'Color','Green');
        hold on
        plot3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4),'LineWidth',2,'Color','Red');
        hold on
        % Impulse at begining
        quiver3(Y_spacecraft_launch_impact(1,2),Y_spacecraft_launch_impact(1,3),Y_spacecraft_launch_impact(1,4), V_impulse(1), V_impulse(2), V_impulse(3), 'b', 'AutoScaleFactor', 0.1);
        % Sun
        scatter3(0,0,0,'red')

        % Create MP4
        v = VideoWriter([pathTosaveAnimation,num2str(PHA),'_',save_type,'_transfer_animation','.mp4'],'MPEG-4');


        v.FrameRate = 30;
        open(v)

        for i = 1 : length(Y_spacecraft_launch_impact)

            if exist('Earth_point', 'var') && ~isempty(Earth_point)
                delete(Earth_point);
                delete(PHA_point);
                delete(Spacecraft_point);
                delete(PHA_orbit);
                delete(Earth_orbit);
                delete(Spacecraft_orbit)
            end

            % Point
            Earth_point = scatter3(Y_Earth_launch_impact(1,i), Y_Earth_launch_impact(2,i),Y_Earth_launch_impact(3,i), 100, 'b', 'filled');
            PHA_point = scatter3(Y_PHA_launch_impact(1,i), Y_PHA_launch_impact(2,i),Y_PHA_launch_impact(3,i), 100, 'g', 'filled');
            Spacecraft_point = scatter3(Y_spacecraft_launch_impact(i,2), Y_spacecraft_launch_impact(i,3), Y_spacecraft_launch_impact(i,4), 100, 'r', 'filled');

            % Orbit
            Earth_orbit = plot3(Y_Earth_launch_impact(1,1:i),Y_Earth_launch_impact(2,1:i),Y_Earth_launch_impact(3,1:i),'LineWidth',2,'Color','Blue');
            hold on
            PHA_orbit = plot3(Y_PHA_launch_impact(1,1:i),Y_PHA_launch_impact(2,1:i),Y_PHA_launch_impact(3,1:i),'LineWidth',2,'Color','Green');
            hold on
            Spacecraft_orbit = plot3(Y_spacecraft_launch_impact(1:i,2),Y_spacecraft_launch_impact(1:i,3),Y_spacecraft_launch_impact(1:i,4),'LineWidth',2,'Color','Red');
            hold on
            % Accerelation
            Accerelation = quiver3(Y_spacecraft_launch_impact(i,2), Y_spacecraft_launch_impact(i,3), Y_spacecraft_launch_impact(i,4), a_low_thrust(1,i)', a_low_thrust(2,i)', a_low_thrust(3,i)', 'r', 'AutoScaleFactor', 0.3);
            lg = legend('Earth','PHA','Spacecraft Trajectory');
            lg.Units = 'inches';
            lg.Position = [12.5,8.5,0,0];

            xlabel('X');
            xlim([-3 3]);
            ylabel('Y');
            ylim([-2 2]);
            zlabel('Z');
            zlim([-3 1.2]);
            title(['Impulse-Low-Thrust acceleration to deflect ',num2str(PHA),PHA_name{PHA},' with ',save_type, ' strategy']);
            grid on;
            view(-40,70);

            % get the frame to video
            frame = getframe(gcf);
            if i == 1
                % Write 10 times to slow the animation down
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
            else
                % Write 3 times to slow the animation down
                writeVideo(v,frame)
                writeVideo(v,frame)
                writeVideo(v,frame)
            end


        end
        close(v)
        close
% -------------------------------------------------------------------------

    end
end





