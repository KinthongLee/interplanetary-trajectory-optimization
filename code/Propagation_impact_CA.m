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
% This function will propagate the PHA state from date of impact to date 
% of closest-approach (CA), and calculate the deflection distance. The
% details of calculation are in paper:
% Input:
% 1. delta_v to PHA after impact, 2. Particles' Parameters, 3. problem
% parameters
% Output:
% 1. DeflectionDistance
%--------------------------------------------------------------------------
function DeflectionDistance = Propagation_impact_CA(delta_v,Particle_Params,problem)
% Impact date
[years0,months0, days0, fds0] = iauJd2cal( 2400000.5, Particle_Params.impact_date);
[hours0, minutes0, secs0] = fd_to_hms(fds0);
[years0, months0, days0, hours0, minutes0, secs0] = fix_seconds(years0, months0, days0, hours0, minutes0, secs0);
secs0 = round(secs0,3);
months0_string = monthToString(months0);
str = sprintf(' %s %g , %g %g:%g:%g', months0_string, days0, years0, hours0, minutes0, secs0);
% Obtain ET (Ephemeris Time )
et_impact = cspice_str2et(str);

% SOI date, the date PHA enter Earth's Sphere of Influence
[years0,months0, days0, fds0] = iauJd2cal( 2400000.5, problem.date_SOI);
[hours0, minutes0, secs0] = fd_to_hms(fds0);
[years0, months0, days0, hours0, minutes0, secs0] = fix_seconds(years0, months0, days0, hours0, minutes0, secs0);
secs0 = round(secs0,3);
months0_string = monthToString(months0);
str = sprintf(' %s %g , %g %g:%g:%g', months0_string, days0, years0, hours0, minutes0, secs0);
% Obtain ET (Ephemeris Time )
et_SOI = cspice_str2et(str);

% Obtain PHA state at moment of impact
[Y, ~] = cspice_spkezr(problem.Target, et_impact, 'J2000', 'NONE', 'SUN');
Y = Y.*1000; % m/s

% Obtain Earth state at the moment of PHA enter Earth's SOI
[Y_EARTH_SOI_heliocentric, ~] = cspice_spkezr('EARTH', et_SOI, 'J2000', 'NONE', 'SUN');
Y_EARTH_SOI_heliocentric = Y_EARTH_SOI_heliocentric.*1000;

% Add delta_V to state of PHA at the moment of impact
Y(4:6,1) = Y(4:6) + delta_v;


% ----------------------- Propagation--------------------------------------
% Impact to SOI, in Heliocentric
delta_t_impact_SOI = et_SOI - et_impact;
Y_impact_SOI = Propagation_heliocentric(Y, delta_t_impact_SOI, Particle_Params.impact_date);

% SOI to CA + 1 days, in Geocentric
delta_t_SOI_CA =  ( problem.date_CA + 1 - problem.date_SOI ) * 86400;
Y_SOI_geocentric = Y_impact_SOI(end,2:7)' - Y_EARTH_SOI_heliocentric;
Y_SOI_CA_geocentric = Propagation_Geocentric(Y_SOI_geocentric ,  delta_t_SOI_CA , problem.date_SOI);
% -------------------------------------------------------------------------

% Distance between PHA & Earth
height = zeros(length( Y_SOI_CA_geocentric),1);
for i = 1 : length( Y_SOI_CA_geocentric)
    height(i) = norm(Y_SOI_CA_geocentric(i,2:4)');
end

% Deflection Distance = Minimum height of PHA to Earth - PHA's CA distance
DeflectionDistance = min(height) - problem.CA_distance ;
end