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
% This code will plot 2 figures as used in the paper:
% 1. PSO Interception Error History
% 2. PSO Deflection Distance History
%--------------------------------------------------------------------------


clear


%--------------------------------------------------------------------------
%------------------- Read data from PHA_table.xlsx ------------------------
%--------------------------------------------------------------------------
filename = 'PHA_table.xlsx';
data = readtable(filename, 'PreserveVariableNames', true);
PHA_name = string(data.Object);
PHA_bsp = string(data.BSP_file_name);


%--------------------------------------------------------------------------
% ---------------------------  Add Path  ----------------------------------
%--------------------------------------------------------------------------
currentDir = fileparts(mfilename('fullpath'));
addpath(genpath(currentDir))

% different path separators for Window(\), MacOS(/), and Linux(/) 
if ispc
    pathTomatFile = [fullfile(currentDir, 'output_result'),'\','matfile','\'];
    pathTosave_intercept = [fullfile(currentDir, 'output_result','figure','PSO_history','Intercept_error'),'\'];
    pathTosave_deflection = [fullfile(currentDir, 'output_result','figure','PSO_history','Deflection'),'\'];
else
    pathTomatFile = [fullfile(currentDir, 'output_result'),'/','matfile','/'];
    pathTosave_intercept = [fullfile(currentDir, 'output_result','figure','PSO_history','Intercept_error'),'/'];
    pathTosave_deflection = [fullfile(currentDir, 'output_result','figure','PSO_history','Deflection'),'/'];
end




% -------------------------------------------------------------------------
%------------------------- Code starts here -------------------------------
% -------------------------------------------------------------------------
% PHA number, as label in PHA_table.xlsx
for PHA = 1:1
    for is_BIP = 0 : 0
        problem.is_BIP = is_BIP;
        problem.Target_name = char(PHA_name(PHA));
        problem.Target = char(PHA_bsp(PHA)); % Target, 2099942 represents Apophis
        if problem.is_BIP
            FullMATPathName = [pathTomatFile,problem.Target_name, '_BIP.mat'];
            save_type = 'BIP';
        else
            FullMATPathName = [pathTomatFile,problem.Target_name, '_COG.mat'];
            save_type = 'COG';
        end

        % Check is .mat file available
        if isfile(FullMATPathName)
            % Load the mat file
            load(FullMATPathName);
            disp('Mat file read successfully....');
        else
            disp('simulation result .mat file not found, please run main.m code to generate it')
        end

        % Store the data from mat file
        length_i = size(out_final.pop,1);

        % iteration might stop early, so need to calculate it
        j= 1;
        while ~isempty(out_final.pop(1,j).Position)
            length_j = j;
            j = j + 1;
            if j > size(out_final.pop,2)
                break
            end
        end
        Interception = zeros(length_i,length_j);
        DeflectionDistance = zeros(length_i,length_j);
        cost = zeros(length_i,length_j);
        X = zeros(1,length_j);
        Y = zeros(1,length_j);
        for j = 1 : length_j
            for i = 1 : length_i
                % Iteration might stop early, if so, the iteration length
                % will be less than length_j, break it
                if isempty(out_final.pop(i,j).Cost)
                    break 
                end
                A = out_final.pop(i,j);
                cost(i,j) = A.Cost;
                Interception(i,j) = A.Interception;
                DeflectionDistance(i,j) = A.DeflectionDistance;
                Position(i,:,j) = A.Position;
            end
            [~,b] = min(cost(:,j));
            X(1,j) = Interception(b,j);
            Y(1,j) = DeflectionDistance(b,j)/1000;
            Z(j,:) = Position(b,:,j);
        end
        
% -------------------------- Start Plotting -------------------------------
        % Interception Error
        fig1 = figure(1);
        scatter(1:length_j,X,'red')
         set(gca, 'YScale', 'log');
         hold on
        
        title(sprintf(['Histogram of Interception Errors Across Iterations \n',num2str(PHA) , PHA_name{PHA}, ' with ', save_type,' strategy']))
        
        ylabel('Interception Error (m)')
        xlabel('Iterations')
        xlim([0 length(X)])
        exportgraphics(fig1, [pathTosave_intercept,num2str(PHA),save_type,'_PSO_Intercept_error_History.png']);
        close
        
        % Deflection Distance 
        fig2 = figure(2);
        scatter(1:length_j,Y,'red')
        set(gca, 'YScale', 'log');
        hold on
        
        title(sprintf(['Histogram of Deflection Distance Across Iterations \n',num2str(PHA) , PHA_name{PHA}, ' with ', save_type,' strategy']))
        ylabel('Deflection Distance (km)')
        xlabel('Iterations')
        xlim([0 length(Y)])
        exportgraphics(fig2, [pathTosave_deflection,num2str(PHA),save_type,'_PSO_DeflectionDistance_History.png']);
        close

    end
end

