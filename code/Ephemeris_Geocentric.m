function Eph = Ephemeris_Geocentric(Y0,  delta_t_SOI_CA, start_date)

options = rdpset('RelTol',1e-12,'AbsTol',1e-14);

% Step size set to 1 seconds
[t,yout] = ode113(@(t,y) Accel_Geocentric(t, y,start_date),(0 : 1 :   delta_t_SOI_CA),Y0,options);
Eph(:,1) = t;
Eph(:,2:7) = yout;