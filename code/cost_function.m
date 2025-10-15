%==========================================================================
% Interplanetary Trajectory Optimization for Planetary Defense Mission
% (Asteroid Kinetic Delfection)
% 
% Kin Thong Lee
% Sept 22 2025
%==========================================================================
% You are free to use and modify the code, but you MUST cite the following
% papers:
%
% Lee, Kinthong, Zhengqing Fang, and Zhaokui Wang. "Investigation of the 
% incremental benefits of eccentric collisions in kinetic deflection of 
% potentially hazardous asteroids." Icarus 425 (2025): 116312.
%
% Feels free to contact me! my email bellow:
% ktlee3819@gmail.com
%==========================================================================
% This is the costfunction (or objective function). The details of
% calculation is described in the paper:
% Given an input of particle's position (1x17 vector), output: 
% 1. Cost Value, 2. Interception error, 3. DeflectionDistance,
% 4. InterceptionAngle, 5. Remaining Mass at Impact Moment
%--------------------------------------------------------------------------
function [cost, Interception, DeflectionDistance, InterceptionAngle, ImpactMass]= cost_function(x,problem)
% Extract value from position vector
Particle_Params.thetaIN_coeff = x(1:6);
Particle_Params.thetaOUT_coeff = x(7:12);
Particle_Params.v_inf_abs = x(13);
Particle_Params.RA = x(14);
Particle_Params.DEC = x(15);
Particle_Params.launch_date = x(16);
Particle_Params.impact_date = x(17);

% Units
AU = 149597870691; % m/AU
TU = 5022548.032; % s/TU


% Launch Date
[years_launch,months_launch, days_launch, fds_launch] = iauJd2cal( 2400000.5, Particle_Params.launch_date);
[hours_launch, minutes_launch, secs_launch] = fd_to_hms(fds_launch);
[years_launch, months_launch, days_launch, hours_launch, minutes_launch, secs_launch] = fix_seconds(years_launch, months_launch, days_launch, hours_launch, minutes_launch, secs_launch);
secs_launch = round(secs_launch,3);
months_launch_string = monthToString(months_launch);
str = sprintf(' %s %g , %g %g:%g:%g', months_launch_string, days_launch, years_launch, hours_launch, minutes_launch, secs_launch);
% Obtain ET (Ephemeris Time )
et_launch = cspice_str2et(str);
% Obtain initial conditions from SPICE
[Y_launch, ~] = cspice_spkezr('EARTH', et_launch, 'J2000', 'NONE', 'SUN');
% State of spacecraft = State of the Earth
r = [Y_launch(1);Y_launch(2);Y_launch(3)];
v = [Y_launch(4);Y_launch(5);Y_launch(6)];

% Before adding impulse to spacecraft, 
% Construct the rotation matrix to transform to heliocentric frame
% Unit vectors of the local orbital frame
r_hat = r / norm(r); % Radial unit vector
h = cross(r, v);     % Specific angular momentum vector
h_hat = h / norm(h); % Normal unit vector (orbital plane normal)
t_hat = cross(h_hat, r_hat); % Tangential unit vector

% Rotation matrix (orbital frame to heliocentric frame)
R = [r_hat, t_hat, h_hat]; % Each column is a unit vector

% Transform impulse to heliocentric frame
delta_v = R * [cos(Particle_Params.DEC)*cos(Particle_Params.RA);cos(Particle_Params.DEC)*sin(Particle_Params.RA);sin(Particle_Params.DEC)]; % Transform orbital to heliocentric frame

% Add impulse (v_inf) to initial velocity
V_impulse = delta_v .* Particle_Params.v_inf_abs;
Y_launch(4:6) = Y_launch(4:6) + V_impulse;
% Convert to m, m/s
Y_launch = Y_launch.*1000; % m, m/s
% Convert to AU, AU/s
Y_launch = Y_launch ./ AU; % AU, AU/s
% Convert velocity to AU/TU
Y_launch(4:6) = Y_launch(4:6) .* TU; % AU/TU


% Impact Date
[years_impact,months_impact, days_impact, fds_impact] = iauJd2cal( 2400000.5, Particle_Params.impact_date);
[hours_impact, minutes_impact, secs_impact] = fd_to_hms(fds_impact);
[years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact] = fix_seconds(years_impact, months_impact, days_impact, hours_impact, minutes_impact, secs_impact);
secs_impact = round(secs_impact,3);
months_impact_string = monthToString(months_impact);
str = sprintf(' %s %g , %g %g:%g:%g', months_impact_string, days_impact, years_impact, hours_impact, minutes_impact, secs_impact);
% Obtain ET (Ephemeris Time )
et_impact = cspice_str2et(str);
% Obtain PHA's conditions at impact date from SPICE
[Y_impact, ~] = cspice_spkezr(problem.Target, et_impact, 'J2000', 'NONE', 'SUN');
Y_impact = Y_impact.*1000; % m, m/s
Y_impact = Y_impact ./ AU; % AU, AU/s
Y_impact(4:6) = Y_impact(4:6) .* TU; % AU/TU



% Propagation from launch date to impact date, with impulse and low-thrust 
% maneuver. To Obtain the v_spacecraft and interception error at impact moment

% Make Sure Impact date is larger than Launch Date (Even we did handle this problem 
% at PSO function, but during the PSO iteration calculation, the impact date
% still could be smaller than launch date, so need a if decision to treat this 
% situation.)
if et_impact > et_launch
% ---------------- Propagation from launch date to impact date ------------
    tspan = [0, (et_impact - et_launch)/TU];
    % Propagate spacecraft from launch date to impact date
    Y_spacecraft_launch_impact = Propagation_lowthrust_heliocentric(Y_launch, tspan, Particle_Params,problem);
    
    % Calculate interception error
    % Interception error = || r_spacecraft - r_PHA ||
    constrain_for_collision = norm(Y_spacecraft_launch_impact(end,2:4)' - Y_impact(1:3));
    % Convert to meter
    Interception = constrain_for_collision * AU; % m

    % Calculate Interception Angle
    InterceptionAngle = acosd( dot(Y_spacecraft_launch_impact(end,5:7)' ,  Y_impact(4:6)) /norm(Y_spacecraft_launch_impact(end,5:7)') /norm(Y_impact(4:6)) );

    % Calculate remaining mass
    % After impulse
    mass_after_depart_Earth = problem.m0 * exp( - Particle_Params.v_inf_abs / problem.c_impulse ); % kg
    % After low-thrust
    mass_at_impact = mass_after_depart_Earth * exp( - problem.n0 / (problem.c_lowthrust*1000) * (et_impact - et_launch) ); % kg
    ImpactMass = mass_at_impact;
% -------------------------------------------------------------------------



% ---------------------------- Impact model -------------------------------
    % Calculate delta_v after impact
    % Vertices of 3D model
    vertices = problem.vertices;
    % Faces of 3D model
    faces = problem.faces;

    % Collision
    if problem.is_BIP 
        % BIP stategy
        [delta_v,~,~] = get_delta_v_from_momentum(Y_spacecraft_launch_impact(end,5:7)', Y_impact(4:6), mass_at_impact ,problem,vertices,faces);
        delta_v_to_PHA_at_impact = delta_v' .* AU ./ TU; % m, m/s, switch from 1x3 to 3x1
    else
        % COG strategy
        [~,delta_v,~] = get_delta_v_from_momentum(Y_spacecraft_launch_impact(end,5:7)', Y_impact(4:6), mass_at_impact ,problem,vertices,faces);
        delta_v_to_PHA_at_impact = delta_v .* AU ./ TU; % m, m/s, NO NEED switch, it is already 3x1
    end
% -------------------------------------------------------------------------
    

% --------- Propagation from impact date to Closest-Approach Date ---------
    % Propagation from impact date to CA date to calculate deflection distance
    % The proppagation in the function is divided into two step:
    % 1. from Impact date to Earth SOI : in heliocentric, a step size of 1 day (86400secs)
    % 2. from Earth SOI to Ca date + 1 day : in geocentric, a step size of 1 seconds
    DeflectionDistance = Propagation_impact_CA(delta_v_to_PHA_at_impact,Particle_Params,problem); % m
% -------------------------------------------------------------------------
    
    % Cost value (Objective value)
    cost = (-DeflectionDistance ) +  100*Interception; % m;

else
    % If impact date < launch date, invalid, make everything Inf
    cost = inf;
    Interception = inf;
    DeflectionDistance = inf;
    InterceptionAngle = inf;
    ImpactMass = inf;
end
    

    
end





