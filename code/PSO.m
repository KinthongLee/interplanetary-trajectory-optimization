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
% This is the code for Particle Swarm Optimization (PSO)
%--------------------------------------------------------------------------
function out = PSO(problem, params)

    % Problem Definiton
    CostFunction = problem.CostFunction;  % Cost Function
    nVar = problem.nVar;        % Number of Unknown (Decision) Variables
    VarSize = [1 nVar];         % Matrix Size of Decision Variables
    PositionMin = problem.PositionMin;	% Lower Bound of Decision Variables
    PositionMax = problem.PositionMax;    % Upper Bound of Decision Variables
    VelocityMin = problem.VelocityMin;
    VelocityMax = problem.VelocityMax;

    % Parameters of PSO
    MaxIt = params.MaxIt;   % Maximum Number of Iterations
    nPop = params.nPop;     % Population Size (Swarm Size)
    ci = params.ci; 
    cc = params.cc;
    cs = params.cs;

    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    

    % Initialization
    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Interception = [];
    empty_particle.DeflectionDistance = [];
    empty_particle.InterceptionAngle = [];
    empty_particle.ImpactMass = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];
    empty_particle.Best.Interception = [];
    empty_particle.Best.DeflectionDistance = [];
    empty_particle.Best.InterceptionAngle = [];
    empty_particle.Best.ImpactMass = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, MaxIt);

    % Create Temporary Stored variable
    % This temporary variable is created due to the limitation of 
    % parfor loops. We cannot perform the following addition:
    % Position(i+1th iterations) = Position(ith ite) + Valocity(ith ite)
    % Parallel loops doesn't allow it, thus I used temporary to store the previous iterations value
    temporary = repmat(empty_particle,nPop,1);

    % Initialize Global Best
    GlobalBest.Cost = inf;
    GlobalBest.Interception = inf;
    GlobalBest.DeflectionDistance = inf;
    GlobalBest.InterceptionAngle = inf;
    GlobalBest.ImpactMass = inf;
    GlobalBest.Position = [];

    % Initialize Population Members
    parfor i=1:nPop
        % Generate Random Position
        particle(i,1).Position = unifrnd(PositionMin, PositionMax, VarSize);
        % If Impact date < Launch date
        if particle(i,1).Position(1,17) <= particle(i,1).Position(1,16)
            % Regenerate a random value between launch date and max value allowed for impact date
            particle(i,1).Position(1,17) = unifrnd( particle(i,1).Position(1,16) , PositionMax(1,17) ); 
        end
        
        % Store the value for temporary
        temporary(i).Position = particle(i,1).Position;
     
        % Initialize Velocity
        particle(i,1).Velocity = zeros(VarSize);
        temporary(i).Velocity = particle(i,1).Velocity;

        % Evaluation, calculate the cost value by using CostFunction (Objective Function)
        [particle(i,1).Cost, particle(i,1).Interception, particle(i,1).DeflectionDistance, particle(i,1).InterceptionAngle, particle(i,1).ImpactMass] = CostFunction(particle(i,1).Position);
        temporary(i).Cost = particle(i,1).Cost;

        % Update the Personal Best
        particle(i,1).Best.Position = particle(i,1).Position;
        particle(i,1).Best.Cost = particle(i,1).Cost;
        particle(i,1).Best.Interception = particle(i,1).Interception;
        particle(i,1).Best.DeflectionDistance = particle(i,1).DeflectionDistance;
        particle(i,1).Best.InterceptionAngle = particle(i,1).InterceptionAngle;
        particle(i,1).Best.ImpactMass = particle(i,1).ImpactMass;
        temporary(i).Best.Position = particle(i,1).Position;
        temporary(i).Best.Cost = particle(i,1).Cost;
        temporary(i).Best.Interception = particle(i,1).Interception;
        temporary(i).Best.DeflectionDistance = particle(i,1).DeflectionDistance;
        temporary(i).Best.InterceptionAngle = particle(i,1).InterceptionAngle;
        temporary(i).Best.ImpactMass = particle(i,1).ImpactMass;
    end

       
    % Manually modify particles to specific value?
    if params.IsSpecificPosition
        for i = 1 : size(params.specificPosition,1)
            % Generate Random Solution
            particle(i,1).Position = params.specificPosition(i,:);
            temporary(i).Position = particle(i,1).Position;
            
            % Initialize Velocity
            particle(i,1).Velocity = zeros(VarSize);
            temporary(i).Velocity = particle(i,1).Velocity;
            % Evaluation
            [particle(i,1).Cost, particle(i,1).Interception, particle(i,1).DeflectionDistance, particle(i,1).InterceptionAngle, particle(i,1).ImpactMass] = CostFunction(particle(i,1).Position);
            temporary(i).Cost = particle(i,1).Cost;


            % Update the Personal Best
            particle(i,1).Best.Position = particle(i,1).Position;
            particle(i,1).Best.Cost = particle(i,1).Cost;
            particle(i,1).Best.Interception = particle(i,1).Interception;
            particle(i,1).Best.DeflectionDistance = particle(i,1).DeflectionDistance;
            particle(i,1).Best.InterceptionAngle = particle(i,1).InterceptionAngle;
            particle(i,1).Best.ImpactMass = particle(i,1).ImpactMass;
            temporary(i).Best.Position = particle(i,1).Position;
            temporary(i).Best.Cost = particle(i,1).Cost;
            temporary(i).Best.Interception = particle(i,1).Interception;
            temporary(i).Best.DeflectionDistance = particle(i,1).DeflectionDistance;
            temporary(i).Best.InterceptionAngle = particle(i,1).InterceptionAngle;
            temporary(i).Best.ImpactMass = particle(i,1).ImpactMass;
        end
        
    end

        % Update Global Best
        for i = 1 : nPop
            if particle(i,1).Best.Cost < GlobalBest.Cost
                 GlobalBest = particle(i,1).Best;
            end
        end




% ------------------------- Main Loop of PSO ------------------------------
    for it=2:MaxIt
        parfor i=1:nPop
            % Update Velocity
            particle(i,it).Velocity = ...
                  ci*(1+rand())/2*temporary(i).Velocity ...
                + cc*rand()*(temporary(i).Best.Position - temporary(i).Position) ...
                + cs*rand()*(GlobalBest.Position - temporary(i).Position);


            % Apply Velocity Limits
            particle(i,it).Velocity = max(particle(i,it).Velocity, VelocityMin);
            particle(i,it).Velocity = min(particle(i,it).Velocity, VelocityMax);
            % Update Temporary value
            temporary(i).Velocity = particle(i,it).Velocity;
            

            % Update Position
            particle(i,it).Position = temporary(i).Position + particle(i,it).Velocity;
            % Apply Lower and Upper Bound Limits
            particle(i,it).Position = max(particle(i,it).Position, PositionMin);
            particle(i,it).Position = min(particle(i,it).Position, PositionMax);
            % Update Temporary value
            temporary(i).Position = particle(i,it).Position;

            % Evaluation, calculate the cost value by using CostFunction (Objective Function)
            [particle(i,it).Cost, particle(i,it).Interception, particle(i,it).DeflectionDistance, particle(i,it).InterceptionAngle, particle(i,it).ImpactMass] = CostFunction(particle(i,it).Position);
            temporary(i).Cost = particle(i,it).Cost;

            % Update Personal Best
            if particle(i,it).Cost < temporary(i).Best.Cost
                particle(i,it).Best.Position = particle(i,it).Position;
                particle(i,it).Best.Cost = particle(i,it).Cost;
                particle(i,it).Best.Interception = particle(i,it).Interception;
                particle(i,it).Best.DeflectionDistance = particle(i,it).DeflectionDistance;
                particle(i,it).Best.InterceptionAngle = particle(i,it).InterceptionAngle;
                particle(i,it).Best.ImpactMass = particle(i,it).ImpactMass;
                temporary(i).Best.Position = particle(i,it).Position;
                temporary(i).Best.Cost = particle(i,it).Cost;
                temporary(i).Best.Interception = particle(i,it).Interception;
                temporary(i).Best.DeflectionDistance = particle(i,it).DeflectionDistance;
                temporary(i).Best.InterceptionAngle = particle(i,it).InterceptionAngle;
                temporary(i).Best.ImpactMass = particle(i,it).ImpactMass;
            end
        end

        % Update Global Best
        for j = 1 : nPop
            if particle(j,it).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(j,it).Best;
            end   
        end

        % Display Iteration Information?
        if ShowIterInfo
            fprintf('Iteration %i: Best Cost = %d, Best Interception = %d (m), Best Deflection Distance = %d (km)\n', it, GlobalBest.Cost, GlobalBest.Interception, GlobalBest.DeflectionDistance/1000)
            disp(vpa(GlobalBest.Position,4))
        end

        % From my experience, if the deflection distance remains negative after more than 80 iterations,
        % then regardless of how many additional iterations are performed, 
        % the computation will not yield a positive deflection distance. 
        % In such cases, the PSO process should be terminated. 
        % Please rerun the calculation; 
        % if repeated runs still fail to produce a positive deflection distance, 
        % it is likely that the asteroid cannot be deflected
        % for example, as demonstrated in my paper with the case of 2017 VW13.
        if it > 80
            if GlobalBest.DeflectionDistance < 0
                disp(it)
                break
            end
        end

        % % If the deflection distance is greater than 0 and the interception error is less than 0.01 m,
        % % the computation is complete and should be terminated.
        % if GlobalBest.DeflectionDistance > 0 && GlobalBest.Interception < 0.01
        %     disp(it)
        %     break
        % end

    end

    % Output
    out.pop = particle;
    out.BestSol = GlobalBest;

    
end