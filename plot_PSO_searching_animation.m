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
% This code will plot two PSO searching animations:
% 1. Trajectory searching
% 2. DeflectionDistance - Intercept Error axis, PSO searching history
%--------------------------------------------------------------------------


clear

%--------------------------------------------------------------------------
% ---------------------------  Add Path  ----------------------------------
%--------------------------------------------------------------------------
currentDir = fileparts(mfilename('fullpath'));
addpath(genpath(currentDir))

% different path separators for Window(\), MacOS(/), and Linux(/) 
if ispc
    pathTomatFile = [fullfile(currentDir, 'output_result','matfile'),'\'];
    pathToPSOSearchMatfile = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_data'),'\'];
    pathTosave_trajectory_searching = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_trajectory_searching'),'\'];
    pathTosave_history_searching = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_history_searching'),'\'];
    
else
    pathTomatFile = [fullfile(currentDir, 'output_result','matfile'),'/'];
    pathToPSOSearchMatfile = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_data'),'/'];
    pathTosave_trajectory_searching = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_trajectory_searching'),'/'];
    pathTosave_history_searching = [fullfile(currentDir, 'output_result','animation','PSO_searching','PSO_history_searching'),'/'];
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



% -------------------------------------------------------------------------
%------------------- Read data from PHA_table.xlsx ------------------------
% -------------------------------------------------------------------------
% Read PHA_table.xlsx
data = readtable('PHA_table');

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
monthFullNames = {'January', 'February', 'March', 'April', 'May', ...
                  'June', 'July', 'August', 'September', 'October', ...
                  'November', 'December'};

% Create containers.Map to map month abbreviations to month numbers 
% and full names
monthNumMap = containers.Map(monthAbbreviations, monthNumbers);
monthNameMap = containers.Map(monthAbbreviations, monthFullNames);

% Convert month abbreviations to numbers and full names
month_CA_num = cell2mat(values(monthNumMap, tokensMatrix(:,2)));
month_CA_string = values(monthNameMap, tokensMatrix(:,2));
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
    for is_BIP = 1 : 1
        % Problem defination
        problem.is_BIP = is_BIP;
        problem.Target_name = char(PHA_name(PHA));
        problem.Target = char(PHA_bsp(PHA)); % Target, 2099942 represents Apophis

        if problem.is_BIP
            FullPSOPathName = [pathToPSOSearchMatfile,problem.Target_name, '_PSOsearching_BIP.mat'];
            FullMATPathName = [pathTomatFile,problem.Target_name, '_BIP.mat'];
            save_type = 'BIP';
        else
            FullPSOPathName = [pathToPSOSearchMatfile,problem.Target_name, '_PSOsearching_COG.mat'];
            FullMATPathName = [pathTomatFile,problem.Target_name, '_COG.mat'];
            save_type = 'COG';
        end

% -------------------- Read / Rerun PSO_searching.mat file ----------------
        % Basically will use orbit propagator run for every single
        % particle's position, and store everything in a .mat file.
        % The .mat file will be read and plot the animation

        % Check is the PSO searching .mat file available
        % If yes, load it
        % If no, rerun and store it
        if isfile(FullPSOPathName)
            % File exists, read it using readtable
            load(FullPSOPathName);
            disp('PSO file read successfully.');
        else
            % File not exist, rerun the process
            disp('PSO_searching file not file, rerun the process....')
                
            % Check is .mat file available
            if isfile(FullMATPathName)
                % Load the mat file
                load(FullMATPathName);
                disp('Mat file read successfully....');

                % Store the data from mat file
                length_i = size(out_final.pop,1);

                % iteration might stop early, so need to calculate it
                j= 1;
                while ~isempty(out_final.pop(1,j).Position)
                    length_ite = j;
                    j = j + 1;
                    if j > size(out_final.pop,2)
                        break
                    end
                end


                Interception = zeros(length_i,length_ite);
                DeflectionDistance = zeros(length_i,length_ite);
                cost = zeros(length_i,length_ite);
                Position = zeros(length_i,17,length_ite);
                BestInterception = zeros(1,length_ite);
                BestDeflectionDistance = zeros(1,length_ite);
                BestPosition = zeros(length_ite,17);

                % Obtain every particle's Interception Error and Deflection
                % Distance
                for ite = 1 : length_ite
                    for i = 1 : length_i
                        A = out_final.pop(i,ite);
                        cost(i,ite) = A.Cost;
                        Interception(i,ite) = A.Interception;
                        DeflectionDistance(i,ite) = A.DeflectionDistance;
                        Position(i,:,ite) = A.Position;
                    end
                    [~,b] = min(cost(:,ite));
                    BestInterception(1,ite) = Interception(b,ite);
                    BestDeflectionDistance(1,ite) = DeflectionDistance(b,ite)/1000;
                    BestPosition(ite,:) = Position(b,:,ite);   
                end

                
                Y_spacecraft_launch_impact = cell(length_i, length_ite); % max_i和max_ite是i和ite的最大值

                % Run Orbit Propagation for every particle's history to
                % obtain their transfer trajectory
                for ite = 1 : length_ite
                    parfor i = 1 : length_i
                        Particle_Params = [];
                        Particle_Params.beta_coeff = Position(i,1:6,ite);
                        Particle_Params.gamma_coeff = Position(i,7:12,ite);
                        Particle_Params.v_inf_abs = Position(i,13,ite);
                        Particle_Params.RA = Position(i,14,ite);
                        Particle_Params.DEC = Position(i,15,ite);
                        Particle_Params.launch_date = Position(i,16,ite);
                        Particle_Params.impact_date = Position(i,17,ite);
        
                        % Convert date from MJday to ET(Ephemeris Time )
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
        
                        % Step 3: Transform acceleration to heliocentric frame
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
                        [Y_impact, ~] = cspice_spkezr(problem.Target, et_impact, 'J2000', 'NONE', 'SUN');
                        Y_impact = Y_impact.*1000; % m, m/s
                        Y_impact = Y_impact ./ AU; % AU, AU/s
                        Y_impact(4:6) = Y_impact(4:6) .* TU; % AU/TU
        
        
        
                        % Make Sure Impact date is larger than Launch Date
                        % Calculate low thrust trasnfer tracjectory to obtain v_sc at impact &
                        % Intercept Distance
                        if et_impact - et_launch > 0.1
                            tspan = [0, (et_impact - et_launch)/TU];
                             Y_launch_impact = Ephemeris_lowthrust_heliocentric(Y_launch, tspan,Particle_Params,problem);
                             Y_spacecraft_launch_impact{i, ite} = Y_launch_impact;
                        else
                           Y_spacecraft_launch_impact{i, ite} = NaN;
                        end
        
                    end

                    % Progressbar
                    showProgress(ite, length_ite)

                end


                % Save the file
                save(FullPSOPathName,'Y_spacecraft_launch_impact','Interception','DeflectionDistance','BestDeflectionDistance','BestInterception','BestPosition')
            else
                disp('simulation result .mat file not found, please run main.m code to generate it')
            end
            disp('PSO file generated, starting to plot....')

        end






% -------------------------------------------------------------------------
%-----------------------------Plot orbit Searching-------------------------
% -------------------------------------------------------------------------
        % Find and plot Earth and PHA orbit
        [years_impact,months_impact, days_impact, fds_impact] = iauJd2cal( 2400000.5, date_SOI(PHA)-365*10);
        [hours_impact, minutes_impact, secs_impact] = fd_to_hms(fds_impact);
        [years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact] = fix_seconds(years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact);
        secs_impact = round(secs_impact,3);
        months_impact_string = monthToString(months_impact);
        str = sprintf(' %s %g , %g %g:%g:%g', months_impact_string, days_impact, years_impact, hours_impact, minutes_impact, secs_impact);
        et_1 = cspice_str2et(str);

        % Earth's Orbit
        [Y_Earth, ~] = cspice_spkezr('Earth', et_1:86400:et_1 + 380*86400, 'J2000', 'NONE', 'SUN');
        Y_Earth = Y_Earth.*1000; % m, m/s
        Y_Earth = Y_Earth ./ AU; % AU, AU/s
        Y_Earth(4:6,:) = Y_Earth(4:6,:) .* TU; % AU/TU

        % PHA's Orbit
        [Y_PHA, ~] = cspice_spkezr(problem.Target, et_1:86400:et_1 + 800*86400, 'J2000', 'NONE', 'SUN');
        Y_PHA = Y_PHA.*1000; % m, m/s
        Y_PHA = Y_PHA ./ AU; % AU, AU/s
        Y_PHA(4:6,:) = Y_PHA(4:6,:) .* TU; % AU/TU


        % Start plotting
        fig = figure(1);
        fig.Units = 'inches';
        fig.Position = [0, 1, 15, 11]; % [left, bottom, width, height]
        % Earth's Orbit
        Earth_orbit = plot3(Y_Earth(1,:),Y_Earth(2,:),Y_Earth(3,:),'blue','LineWidth',5);
        hold on
        % PHA's Orbit
        PHA_orbit = plot3(Y_PHA(1,:),Y_PHA(2,:),Y_PHA(3,:),'green','LineWidth',5);
        hold on

        xlabel('X(AU)');
        xmin = min( min(min(Y_PHA(1:3,:))), min(min(Y_Earth(1:3,:))) ) - 0.5;
        xmax = max( max(max(Y_PHA(1:3,:))), max(max(Y_Earth(1:3,:))) ) + 0.5;
        xlim([xmin xmax]);

        ylabel('Y(AU)');
        ymin = min( min(min(Y_PHA(1:3,:))), min(min(Y_Earth(1:3,:))) ) - 0.5;
        ymax = max( max(max(Y_PHA(1:3,:))), max(max(Y_Earth(1:3,:))) ) + 0.5;
        ylim([ymin ymax]);

        zlabel('Z(AU)');
        zmin = min( min(min(Y_PHA(1:3,:))), min(min(Y_Earth(1:3,:))) ) - 0.1;
        zmax = max( max(max(Y_PHA(1:3,:))), max(max(Y_Earth(1:3,:))) ) + 0.1;
        zlim([zmin zmax]);
        grid on;
        view(-40,70);


        % Create MP4
        v = VideoWriter([pathTosave_trajectory_searching,num2str(PHA),PHA_name{PHA}, save_type, '_PSO','_Trajectory_searching','.mp4'],'MPEG-4');
        v.FrameRate = 30;
        open(v)

        for ite = 1 : size(Y_spacecraft_launch_impact,2)
            % Create a handle, store every plot's name in each iteration
            % in preparation to be deleted in next iteration
            transfer_handles = [];
            for i = 1 : size(Y_spacecraft_launch_impact,1)
                if ~isnan(Y_spacecraft_launch_impact{i,ite})
                    % Plot ith particle's trajectory
                    transfer_orbit = plot3(Y_spacecraft_launch_impact{i,ite}(:,2),Y_spacecraft_launch_impact{i,ite}(:,3),Y_spacecraft_launch_impact{i,ite}(:,4),'red','LineWidth',1);
                    hold on
                    transfer_handles = [transfer_handles, transfer_orbit];
                end
            end

            title(sprintf(['Deflecting',num2str(PHA) , PHA_name{PHA}, ' with ', save_type, ...
                ' strategy \n Iteration: %i \n Best Deflection Distance: %dm Best Interception Error: %dm '], ...
                ite, BestDeflectionDistance(1,ite), BestInterception(1,ite)));

            hLeg = legend([Earth_orbit,PHA_orbit, transfer_orbit], 'Earth orbit', 'PHA orbit', 'transfer trajectory');
            hLeg.Position = [0.7558, 0.5662, 0.1007, 0.0445];

            % get the frame to video
            frame = getframe(gcf);
            writeVideo(v,frame)

            if ite ~= size(Y_spacecraft_launch_impact,2)
                % Delete all particles' trajectories, besides the last
                % iteration. Cannot use cla reset because the Earth and PHA
                % Orbit need to stay
                delete(transfer_handles)
            end

        end
        close(v)
        close
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
%----------------------- Plot PSO searching History -----------------------
% -------------------------------------------------------------------------
fig = figure(2);
fig.Units = 'inches';
fig.Position = [0, 1, 15, 11]; % [left, bottom, width, height]


% Create MP4
v = VideoWriter([pathTosave_history_searching,num2str(PHA),PHA_name{PHA}, save_type,'_PSO_search_History','.mp4'],'MPEG-4');
v.FrameRate = 30;
open(v)

xmin = min(Interception(~isinf(Interception(:,:))));
xmax = max(Interception(~isinf(Interception(:,:))));
ymin = min(DeflectionDistance(~isinf(DeflectionDistance(:,:))));
ymax = max(DeflectionDistance(~isinf(DeflectionDistance(:,:))));


for ite = 1:size(Interception,2)
    % Plot
    scatter(Interception(:,ite), DeflectionDistance(:,ite), 15, 'filled','red'); % 或 'r.'
    hold on

    set(gca, 'XScale','log');
    xlabel('Interception Error (m)')
    ylabel('Deflection Distance (m)');
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    grid on;

    title(sprintf(['Deflecting',num2str(PHA) , PHA_name{PHA}, ' with ', save_type, ...
                ' strategy \n Iteration: %i \n Best Deflection Distance: %dm Best Interception Error: %dm '], ...
                ite, BestDeflectionDistance(1,ite), BestInterception(1,ite)));
    

    drawnow;  

    frame = getframe(gcf);
    writeVideo(v, frame);

    hold off


    if ite ~= size(Interception,2)
        % Delete every thing
        cla reset   
    end
end

close(v);
close(fig);
% -------------------------------------------------------------------------
        

    end
end





