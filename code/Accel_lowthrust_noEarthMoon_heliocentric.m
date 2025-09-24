
function dY = Accel_lowthrust_noEarthMoon_heliocentric(t, Y,start_date,Particle_Params,problem)
% Extract value
thetaIN_coeff = Particle_Params.thetaIN_coeff;
thetaOUT_coeff = Particle_Params.thetaOUT_coeff;
c = problem.c_lowthrust;
n0 = problem.n0;
AU = 149597870691; % m/AU
TU = 5022548.032; % s/TU

% Time
MJD_UTC = start_date + t/86400*TU;
[year,month, day, fd] = iauJd2cal( 2400000.5, MJD_UTC);
[hour, minute, sec] = fd_to_hms(fd);
[year, month, day, hour, minute, sec] = fix_seconds(year, month, day, hour, minute, sec);
sec = round(sec,3);
month = monthToString(month);
% Obtain ET (Ephemeris Time )
str = sprintf(' %s %g , %g %g:%g:%g', month, day, year, hour, minute, sec);
et = cspice_str2et(str);


% ---------------------------- Perturbations ------------------------------
% Obtain position vector of planets
% [r_Mercury, ~] = cspice_spkezr('Mercury', et, 'J2000', 'NONE', 'SUN');
% r_Mercury = r_Mercury(1:3)./(AU/1000); % m/s
[r_Venus, ~] = cspice_spkezr('Venus', et, 'J2000', 'NONE', 'SUN');
r_Venus = r_Venus(1:3)./(AU/1000); % AU
[r_Mars, ~] = cspice_spkezr('Mars Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Mars = r_Mars(1:3)./(AU/1000); % AU
[r_Jupiter, ~] = cspice_spkezr('Jupiter Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Jupiter = r_Jupiter(1:3)./(AU/1000); % AU
% [r_Saturn, ~] = cspice_spkezr('Saturn Barycenter', et, 'J2000', 'NONE', 'SUN');
% r_Saturn = r_Saturn(1:3)./(AU/1000); % m/s
% [r_Uranus, ~] = cspice_spkezr('Uranus Barycenter', et, 'J2000', 'NONE', 'SUN');
% r_Uranus = r_Uranus(1:3)./(AU/1000); % m/s
% [r_Neptune, ~] = cspice_spkezr('Neptune Barycenter', et, 'J2000', 'NONE', 'SUN');
% r_Neptune = r_Neptune(1:3)./(AU/1000); % m/s
% [r_Pluto, ~] = cspice_spkezr('Pluto Barycenter', et, 'J2000', 'NONE', 'SUN');
% r_Pluto = r_Pluto(1:3)./(AU/1000); % m/s


% Gravitational Coefficient, in AU^3/TU^2 units
mu_Sun = 1;
% mu_Earth = 398600441800000/ AU^3 * TU^2;
% mu_Moon = 4902800982000/ AU^3 * TU^2;
mu_Venus = 324858592000000/ AU^3 * TU^2;
mu_Mars = 42828375816000/ AU^3 * TU^2;
mu_Jupiter = 126712764100000000/ AU^3 * TU^2;  


% Sun two body Gravity
a = Accel_two_body(Y(1:3),mu_Sun);

%Planetary perturbations
% a = a + AccelPointMass(Y(1:3),r_Mercury(1:3),22031868551000/ AU^3 * TU^2);
a = a + AccelPointMass(Y(1:3),r_Venus(1:3),mu_Venus);
a = a + AccelPointMass(Y(1:3),r_Mars(1:3),mu_Mars);
a = a + AccelPointMass(Y(1:3),r_Jupiter(1:3),mu_Jupiter);
% a = a + AccelPointMass(Y(1:3),r_Saturn(1:3),37940584841800000/ AU^3 * TU^2);
% a = a + AccelPointMass(Y(1:3),r_Uranus(1:3),5794556400000000/ AU^3 * TU^2);    
% a = a + AccelPointMass(Y(1:3),r_Neptune(1:3),6836527100580000/ AU^3 * TU^2);
% a = a + AccelPointMass(Y(1:3),r_Pluto(1:3),975500000000/ AU^3 * TU^2);
    
% % Non-gravitational forces
% % Solar rad
a = a + AccelsRad(Y(1:3).*AU,149597870700) ./ AU .* TU^2;
%--------------------------------------------------------------------------



% ----------------------- Low-Thrust maneuvers ----------------------------
r = [Y(1);Y(2);Y(3)];
v = [Y(4);Y(5);Y(6)];

%  Define thetaIN & thetaOUT as a function of time
%  make sure thetaIN & thetaOUT coeff is given in row vector
if size(thetaIN_coeff, 1) > size(thetaIN_coeff, 2) % check if it is column vector (Nx1)
    thetaIN_coeff = reshape(thetaIN_coeff, 1, []);  % Reshapes to 1xN
end

if size(thetaOUT_coeff, 1) > size(thetaOUT_coeff, 2) % check if it is column vector (Nx1)
    thetaOUT_coeff = reshape(thetaOUT_coeff, 1, []);  % Reshapes to 1xN
end

for i = 1 : length(thetaIN_coeff)
   thetaIN = thetaIN_coeff(i) * t^(i-1);
   thetaOUT = thetaOUT_coeff(i) * t^(i-1);
end


% Thrust accerelation norm:
a_norm = c * n0 * 1000 / (c * 1000  - n0 * t);

% Thrust components in spacecraft orbit coordinate
a_r = a_norm * cos(thetaOUT) * cos(thetaIN); % Radial acceleration
a_t = a_norm * cos(thetaOUT) * sin(thetaIN); % Tangential acceleration
a_n = a_norm * sin(thetaOUT);             % Normal acceleration
% Acceleration vector in the local orbital frame
a_orbital = [a_r; a_t; a_n];

% Construct the rotation matrix to transform to heliocentric frame
% Unit vectors of the local orbital frame
r_hat = r / norm(r); % Radial unit vector
h = cross(r, v);     % Specific angular momentum vector
h_hat = h / norm(h); % Normal unit vector (orbital plane normal)
t_hat = cross(h_hat, r_hat); % Tangential unit vector

% Rotation matrix (orbital frame to heliocentric frame)
R = [r_hat, t_hat, h_hat]; % Each column is a unit vector

% Transform acceleration to heliocentric frame
a_low_thrust = R * a_orbital; % Transform orbital to heliocentric frame

a = a + a_low_thrust;
%--------------------------------------------------------------------------

dY = [Y(4:6);a];

