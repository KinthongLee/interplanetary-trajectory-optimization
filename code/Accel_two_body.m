function a = Accel_two_body(r, GM)

% Acceleration 
a = -GM * ( r/(norm(r)^3) );
end