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
% Deflecting Potentially Hazardous Asteroid (PHA) by using 
% Impulsive-Continous-Low_Thrust Acceleration maneuver.
% Solve the opmiziming problem using Particle Swarm Optimization (PSO)
% The code is easily to be modified to solve any interplanetary trajectory
% optimization problem, by modifying the costfunction
%--------------------------------------------------------------------------


clear 

% -------------------------------------------------------------------------
% ---------------------------  Add Path  ----------------------------------
% -------------------------------------------------------------------------
currentDir = fileparts(mfilename('fullpath'));
% Add all directory to path
addpath(genpath(currentDir))

% different path separators for Window(\), MacOS(/), and Linux(/) 
if ispc
    pathTosave = [fullfile(currentDir, 'output_result'),'\'];
else
    pathTosave = [fullfile(currentDir, 'output_result'),'/'];
end

% Add path to kernel (to be used for SPICE) 
% Specify the directory where your .bsp files are located
pathTokernel = fullfile(currentDir, 'code', 'kernel');
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
%-----------------------------Load Kernels---------------------------------
% -------------------------------------------------------------------------
% Load the kernels for SPICE 

% Start the SPMD block for parallel execution (parloop)
spmd
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
end


% Load Kernels for non parfor loop
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

% Duration of Warning time
warning_years = 10; % years

% Parameters for collision model
PHA_mass = 2.7e10; % kg
beta = 3.61; % momentum enchancement coefficient
PHA_3D_model = 'Apophis_Model.obj'; % 3D Model of PHA

% Parameters for rocket
c_impulse = 3.09; % Exhaust velocity km/s, I = 315s
c_lowthrust = 29.78  ; % Exhaust velocity km/s I = 3035s
m0 = 10000; % Initial mass, kg
n0 = 18e-6 * 9.81 ; % thrust-to-mass ratio at t0, m/s^2

% Parameters for PSO
PSO_velocity_scale_factor = 0.5; % larger value for a highly dynamic systems.
ci = 1; % inertial weights, it will times (1+rand())/2 in PSO.m
cc = 1.49445; % cognitive weights, i will times rand() in PSO.m
cs = 1.49445; % social (stochastic) weights, it will times rand() in PSO.m





% -------------------------------------------------------------------------
% ----------------------------- Simulation --------------------------------
% -------------------------------------------------------------------------
% Simulation starts here:

% PHA number, as label in PHA_table.xlsx
for PHA = 1 : 1

    % Is SOI date (the date PHA enter Earth's Sphere of Influence) available?
    % if not, calculate it
    if isnan(data.Date_SOI(PHA))
        data.Date_SOI(PHA)  = calculate_SOI(year_CA(PHA),month_CA_num(PHA),date_CA(PHA),PHA_bsp(PHA));
        disp('done')
    end
    

    % COG or BIP strategy? 
    % 0 : 0, only COG strategy
    % 1 : 1, only BIP strategy
    % 0 : 1, for both COG & BIP, calculate COG first, then BIP
    for is_BIP =  0 : 0

% ------------------------ Problem Definiton ------------------------------
        problem.nVar = 17;       % Number of Unknown (Decision) Variables
        % Each particle's position has 17 dimensions, which in order:
        % 5th polynomial coeffients of thetaIN、thetaOUT（12), 
        % Earth departure magnitude, RA angle, and DEC angle,
        % launch date, impact date
        % [thetaIN0, ..., thetaIN5, thetaOUT0, ..., thetaOUT5, Delta_V, RA,
        %  DEC, Mjday_launch, Mjday_impact]

        % launch and impact date
        launch_min = Mjday(year_CA(PHA) - warning_years, 1, 1); 
        launch_max = Mjday(year_CA(PHA), 1, 1); 
        impact_min = launch_min;
        impact_max = launch_max; 
    
        % Calculate the Guess for thetaIN and thetaOUT coefficients
        % thetaIN = thetaIN0 + thetaIN1 * t + thetaIN2* t^2 + ... 
        %           + thetaIN5 * t^5
        % Where the t is in unit of T_TU = 5022548.032 seconds, 
        % 2*PI T_TU = 365.25 days to avoid calculation exceed 
        % matlab's floating-point precision
        thetaIN_coeff_value = [1e2, 1e1, 1e-2, 1e-3, 1e-4, 1e-5];
        thetaOUT_coeff_value = [1e2, 1e1, 1e-2, 1e-3, 1e-4, 1e-5];

        % The elements of position need to be bounded 
        % Min and Max value for Positions
        problem.PositionMin =  [-thetaIN_coeff_value, -thetaOUT_coeff_value, 0, -pi, -pi/2, launch_min, impact_min];  % Lower Bound of Decision Variables
        problem.PositionMax =  [thetaIN_coeff_value, thetaOUT_coeff_value, 20,  pi, pi/2, launch_max, impact_max];   % Upper Bound of Decision Variables
        
        % Min and Max value for Velocity
        problem.VelocityMin = [-thetaIN_coeff_value.*PSO_velocity_scale_factor -thetaOUT_coeff_value.*PSO_velocity_scale_factor -5 -pi/2*PSO_velocity_scale_factor -pi/2*PSO_velocity_scale_factor -365 -365];
        problem.VelocityMax = [thetaIN_coeff_value.*PSO_velocity_scale_factor thetaOUT_coeff_value.*PSO_velocity_scale_factor 5 pi/2*PSO_velocity_scale_factor pi/2*PSO_velocity_scale_factor 365 365];
        
        % Target
        problem.Target_name = char(PHA_name(PHA));
        problem.Target = char(PHA_bsp(PHA)); % etc: 2099942 = Apophis
        problem.mass_PHA = PHA_mass; 
        problem.date_SOI = data.Date_SOI(PHA); % Mjday
        problem.date_CA = Mjday(year_CA(PHA),month_CA_num(PHA),date_CA(PHA));
        problem.CA_distance = data.CADistanceMinimum_km_(PHA) * 1000 ; % m
        
        % Parameters of collision model
        problem.beta = beta;
        problem.is_BIP = is_BIP; 
        problem.PHA_3D_model = PHA_3D_model;
        % Read 3D models
        [problem.vertices, problem.faces] = readObj(problem.PHA_3D_model);

        % Parameter of engine
        problem.c_impulse = c_impulse; 
        problem.c_lowthrust = c_lowthrust  ; 
        problem.m0 = m0; 
        problem.n0 = n0 ; 

        % Cost Function J (or Objective Function)
        % this cost_function.m is store in directory "code"
        problem.CostFunction = @(x) cost_function(x,problem);  

        % Which strategy
        if problem.is_BIP
            strategy_name = 'BIP';
        else
            strategy_name = 'COG';
        end
        
        % Parameters of PSO
        params.ci = ci; 
        params.cc = cc; 
        params.cs = cs; 
        % Parameters for First loop of PSO 
        params.MaxIt = 30;        % Maximum Number of Iterations
        params.nPop = 2;           % Population Size (Swarm Size)
        params.total_loops = 200;      % Total loops 
        params.ShowIterInfo = false; % Flag for Showing Iteration Information
        % Manually give value to one particle? 
        params.IsSpecificPosition = false; 
        %  If true, modify the value below
        % params.specificPosition = [1e+03,-1.07,-4.2e-03,1.6e-04,6.02e-05,-8.5e-06,-2.09e+02,-2.6,-3.8e-03,9.6e-04,6.2e-05,-9.9e-06,1.4,1e+00,-7e-01,5.9193023e+04,5.92893911e+04];
% -------------------------------------------------------------------------


        % Display message
        fprintf('Calculation for PHA %s with %s strategy starts:\n',PHA_name(PHA), strategy_name)
        
        % Create empty
        empty.pop = [];
        empty.BestSol = [];
    
% % ------------------ First loop of PSO, brief searching -------------------
%         fprintf('Calculating first loop for %s: \n', problem.Target_name)
%         out_first = repmat(empty, params.total_loops, 1);
%         DeflectionDistance_first = zeros(params.total_loops,1);
%         % Start timing
%         tic
%         for j = 1 : params.total_loops
%             % PSO run
%             out_first(j) = PSO(problem, params);
%             DeflectionDistance_first(j) = out_first(j).BestSol.DeflectionDistance;
% 
%             % Progressbar
%             showProgress(j, params.total_loops)
%         end
% 
%         % Index for Best Deflection Distance
%         [~,B1] = max(DeflectionDistance_first);
% 
%         % Display Result
%         fprintf('Deflection Distance: %d , Interception Error: %d \n', out_first(B1).BestSol.DeflectionDistance, out_first(B1).BestSol.Interception);
% % -------------------------------------------------------------------------     

        % Specific Modification
            params.IsSpecificPosition = true;
        if is_BIP
            params.specificPosition = Position_BIP(PHA,:); % Value from PHA_table.xlsx
            % params.specificPosition = Position_COG(PHA,:);
        else
            % params.specificPosition = Position_BIP(PHA,:);
            params.specificPosition = [1e+03,-1.07,-4.2e-03,1.6e-04,6.02e-05,-8.5e-06,-2.09e+02,-2.6,-3.8e-03,9.6e-04,6.2e-05,-9.9e-06,1.4,1e+00,-7e-01,5.9193023e+04,5.92893911e+04];
            % params.specificPosition = out_final(1).BestSol.Position;
            % params.specificPosition = out_first(B1).BestSol.Position;
%         params.specificPosition = [-9.9999614705317778629023450776003, -0.49232500328425038427226922976843, 0.0061452204636763749223682395950163, -0.00019420841381091659850338659865088, -0.000022294748096374127836576117811518, -0.0000023349223154250462478271135879915, -9.9999978334439951765943987993523, 0.9999998666087578369499055952474, -0.0076666070152362798817424760500217, 0.00099999950637630558372692668456239, 0.000083219552406479490713059332662738, -0.000003115374897908437071438044960181, 11.222610398950207510893051221501, -1.4750815820997953409232650301419, 0.85699721923916016841360487887869, 54530.390100913362402934581041336, 55691.341553505175397731363773346];
        end
    



% -------------- Final loop of PSO, much detailed searching ---------------
        % Plus minus 30 days of the best Launch Date as new search Date
        % In this final loop, no need so much loops, usually one is enough
        % higher nPop and MaxIt is set for final loop
        
        % % Manually give value to one particle? 
        % params.IsSpecificPosition = true;
        % % Manually add the best position of first loop to one particle
        % params.specificPosition = out_first(B1).BestSol.Position;
        params.total_loops = 1;      % Total loops 
        params.MaxIt = 300;        % Maximum Number of Iterations
        params.nPop = 200;
        params.ShowIterInfo = false; % Flag for Showing Iteration Informatin
        fprintf('Calculation for final loop of PSO for %s starts: \n', problem.Target_name);
        
        out_final = repmat(empty, params.total_loops, 1);
        DeflectionDistance_final = zeros(params.total_loops,1);
        for j = 1 : params.total_loops
            out_final(j) = PSO(problem, params);
            DeflectionDistance_final(j) = out_final(j).BestSol.DeflectionDistance;
        end
        
        % Index for best Deflection Distance
        [~,B2] = max(DeflectionDistance_final);
% -------------------------------------------------------------------------        
    
    

        % Display result
        fprintf('Deflection Distance: %d , Interception Error: %d \n', out_final(B2).BestSol.DeflectionDistance, out_final(B2).BestSol.Interception);
        fprintf('\n')
    
        % Store result into data
        if problem.is_BIP
            filename = [pathTosave,problem.Target_name, '_BIP.mat'];  % Append .mat extension
            % Convert Best Position to cell data
            vectorAsString1 = sprintf('%.30e,', out_final(B2).BestSol.Position(1:16));
            vectorAsString = ['[',vectorAsString1 sprintf('%.30e' , out_final(B2).BestSol.Position(17)),']'];
            cellData = {vectorAsString};  
            data.BestPosition_BIP(PHA) = cellData;
            data.BestDeflectionDistance_BIP(PHA) = out_final(B2).BestSol.DeflectionDistance;
            data.Interception_Error_BIP(PHA) = out_final(B2).BestSol.Interception;
            data.Interception_Angle_BIP(PHA) = out_final(B2).BestSol.InterceptionAngle;
            data.Impact_Mass_BIP(PHA) = out_final(B2).BestSol.ImpactMass;
        else
            filename = [pathTosave,problem.Target_name, '_COG.mat'];  % Append .mat extension
            % Convert Best Position to cell data
            vectorAsString1 = sprintf('%.30e,', out_final(B2).BestSol.Position(1:16));
            vectorAsString = ['[',vectorAsString1 sprintf('%.30e' , out_final(B2).BestSol.Position(17)),']'];
            cellData = {vectorAsString};       
            data.BestPosition_COG(PHA) =cellData;
            data.BestDeflectionDistance_COG(PHA) = out_final(B2).BestSol.DeflectionDistance;
            data.Interception_Error_COG(PHA) = out_final(B2).BestSol.Interception;
            data.Interception_Angle_COG(PHA) = out_final(B2).BestSol.InterceptionAngle;
            data.Impact_Mass_COG(PHA) = out_final(B2).BestSol.ImpactMass;
        end
    
    
        % save.mat file
        save(filename,'out_final','out_first');


        % Write table to .xlsx file (without relative Gain rate)
        writetable(data,'temporary_result.xlsx')
    end
    toc

    % Calculate Relative Gain rate
    data.Gain_percent(PHA) = (data.BestDeflectionDistance_BIP(PHA) - data.BestDeflectionDistance_COG(PHA) ) / data.BestDeflectionDistance_COG(PHA) * 100;
    
    % Write table to .xlsx file (with relative Gain rate)
    writetable(data,'temporary_result.xlsx')
end
    

        
