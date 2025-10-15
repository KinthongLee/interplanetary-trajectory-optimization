function Eph = Ephemeris_lowthrust_heliocentric(Y0, tspan,Particle_Params,problem)

options = rdpset('RelTol',1e-12,'AbsTol',1e-14);

[t,yout] = ode113(@(t,y) Accel_lowthrust_noEarthMoon_heliocentric(t,y,Particle_Params.launch_date,Particle_Params,problem),tspan,Y0,options);
Eph(:,1) = t;
Eph(:,2:7) = yout;

