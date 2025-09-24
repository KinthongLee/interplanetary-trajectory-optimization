
function dY = Accel_Geocentric(t, Y,start_date)
% Iteration time
MJD_UTC = start_date + t/86400;
[year,month, day, fd] = iauJd2cal( 2400000.5, MJD_UTC);
[hour, minute, sec] = fd_to_hms(fd);
[year, month, day, hour, minute, sec] = fix_seconds(year, month, day, hour, minute, sec);
sec = round(sec,3);
month = monthToString(month);
str = sprintf(' %s %g , %g %g:%g:%g', month, day, year, hour, minute, sec);
% Obtain ET (Ephemeris Time )
et = cspice_str2et(str);

% Obtain Sun and Moon positions in Geocentric frame
[r_Sun, ~] = cspice_spkezr('Sun', et, 'J2000', 'NONE', 'Earth');
r_Sun = r_Sun(1:3).*1000; % m/s
[r_Moon, ~] = cspice_spkezr('Moon', et, 'J2000', 'NONE', 'Earth');
r_Moon = r_Moon(1:3).*1000; % m/s


% Gravitational Coefficient
mu_Sun = 132712440041279422464; % m^3/s^2
mu_Earth = 398600441500000;
mu_Moon = 4902800192171;


% Earth two body Gravity
a = Accel_two_body(Y(1:3),mu_Earth);
% J2
a = a + (AccelJ2( Y(1:3), mu_Earth ));
% Sun Perturbation
a = a + AccelPointMass(Y(1:3),r_Sun(1:3),mu_Sun);
% Moon Perturbation
a = a + AccelPointMass(Y(1:3),r_Moon(1:3),mu_Moon);


dY = [Y(4:6);a];