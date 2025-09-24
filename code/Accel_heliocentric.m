function dY = Accel_heliocentric(t, Y, Mjd0_UTC)
% Time
MJD_UTC = Mjd0_UTC+t/86400;
[year,month, day, fd] = iauJd2cal( 2400000.5, MJD_UTC);
[hour, minute, sec] = fd_to_hms(fd);
[year, month, day, hour, minute, sec] = fix_seconds(year, month, day, hour, minute, sec);
sec = round(sec,3);
month = monthToString(month);
% Obtain ET (Ephemeris Time )
str = sprintf(' %s %g , %g %g:%g:%g', month, day, year, hour, minute, sec);
et = cspice_str2et(str);


% Obtain position vector of planets
[r_Moon, ~] = cspice_spkezr('Moon', et, 'J2000', 'NONE', 'SUN');
r_Moon = r_Moon(1:3).*1000; % m/s
[r_Mercury, ~] = cspice_spkezr('Mercury', et, 'J2000', 'NONE', 'SUN');
r_Mercury = r_Mercury(1:3).*1000; % m/s
[r_Venus, ~] = cspice_spkezr('Venus', et, 'J2000', 'NONE', 'SUN');
r_Venus = r_Venus(1:3).*1000; % m/s
[r_Mars, ~] = cspice_spkezr('Mars Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Mars = r_Mars(1:3).*1000; % m/s
[r_Jupiter, ~] = cspice_spkezr('Jupiter Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Jupiter = r_Jupiter(1:3).*1000; % m/s
[r_Saturn, ~] = cspice_spkezr('Saturn Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Saturn = r_Saturn(1:3).*1000; % m/s
[r_Uranus, ~] = cspice_spkezr('Uranus Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Uranus = r_Uranus(1:3).*1000; % m/s
[r_Neptune, ~] = cspice_spkezr('Neptune Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Neptune = r_Neptune(1:3).*1000; % m/s
[r_Pluto, ~] = cspice_spkezr('Pluto Barycenter', et, 'J2000', 'NONE', 'SUN');
r_Pluto = r_Pluto(1:3).*1000; % m/s
[r_Earth, ~] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
r_Earth = r_Earth(1:3).*1000; % m/s


% Sun two body Gravity
a = Accel_two_body(Y(1:3),132712440041279422464);


% Earth Gravity perturbations
    a = a + AccelPointMass(Y(1:3),r_Earth,398600441500000);

% Moon Gravity perturbations
    a = a + AccelPointMass(Y(1:3),r_Moon,4902800192171.3935546875);
% Planetary Gravity perturbations

    a = a + AccelPointMass(Y(1:3),r_Mercury,22031868551000);
    a = a + AccelPointMass(Y(1:3),r_Venus,324858592000000);
    a = a + AccelPointMass(Y(1:3),r_Mars,42828375816000);
    a = a + AccelPointMass(Y(1:3),r_Jupiter,126712764100000000);
    a = a + AccelPointMass(Y(1:3),r_Saturn,37940584841800000);
    a = a + AccelPointMass(Y(1:3),r_Uranus,5794556400000000);    
    a = a + AccelPointMass(Y(1:3),r_Neptune,6836527100580000);
    a = a + AccelPointMass(Y(1:3),r_Pluto,975500000000);

% Non-gravitational forces
% Solar rad
    a = a + AccelsRad(Y(1:3),149597870700);
   % Yarkovsky
    a = a + AccelYarkovsky(Y(1:3),149597870700,Y(4:6));
    




% Relativistic Effects
    a = a + Relativity(Y(1:3),Y(4:6));




dY = [Y(4:6);a];

