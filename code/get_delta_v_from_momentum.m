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
% This is the function of collision model, consider both centric &
% eccentric collision. The details of the calculation are in paper:
% Input: 
% 1. velocity of spacecraft, 2. velocity of PHA, 
% 3. Remaining mass of spacecraft, 4. problem parameters, 
% 5. vertices of PHA 3D model, 6. faces of PHA 3D model
% Output:
% 1. delta_v generated through BIP strategy, in 1x3 form
% 2. delta_v generated through COG strategy, in 3x1 form
% 3. PHA's attitude
%--------------------------------------------------------------------------
function [best_delta_v,delta_v_COM,angles_output] = get_delta_v_from_momentum(v_imp,v_PHA,mass_rocket,problem,vertices,faces)
% Parameters
m = mass_rocket; %kg
M = problem.mass_PHA; %kg
v_r = v_imp-v_PHA;
v_r_norm = v_r / norm(v_r);
beta = problem.beta;

% ------------------------- PHA attitude ----------------------------------
% Define rotation angles, for PHA attitudes
% angles = [rand()*pi;rand()*pi;rand()*pi];
% Below is an attitude that generate average deflection distance, as
% described in papers.
angles = [1.5928340138728427355374606122496, 1.730231432842384009518355014734, 0.28579620163745189653781153538148];

% Build rotation matrix
Rx = [cos(angles(1)),sin(angles(1)),0;-sin(angles(1)),cos(angles(1)),0;0,0,1];
Ry = [cos(angles(2)),0,sin(angles(2));0,1,0;-sin(angles(2)),0,cos(angles(2))];
Rz = [1,0,0;0,cos(angles(3)),sin(angles(3));0,-sin(angles(3)),cos(angles(3))];
R = Rx*Ry*Rz;

% Apply the rotation matrix to the vertices of the model
rotated_vertices = zeros(length(vertices),1);
for i = 1 : length(vertices)
    rotated_vertices(i,1:3) = (R * vertices(i,1:3)')';
end
% -------------------------------------------------------------------------




% ----------------- BIP strategy (Eccentric Strike) -----------------------
% Compute the normal vectors of the every faces
face_normals = cross(rotated_vertices(faces(:,2),:)-rotated_vertices(faces(:,1),:), ...
                     rotated_vertices(faces(:,3),:)-rotated_vertices(faces(:,1),:));
face_normals = face_normals ./ vecnorm(face_normals, 2, 2);
irradiated_vector =  v_r_norm;

% Find all the surface which irradiated by the velocity vector
% Compute the dot product between the face normals and the light direction
dot_products = dot(face_normals, repmat(irradiated_vector', size(faces, 1), 1), 2);

% Find the indices of the irradiated faces
% These faces are the faces facing the spacecraft, which could be 
% potential impact candidates
irradiated_face_indices = dot_products < 0;


% Calculate the tangential_vector_t & theta_impact of every faces
n = -face_normals(irradiated_face_indices,:);  %normal components of face
k = [0,0,1];
t = zeros(length(n),3);  % tangential components of face
theta_impact = zeros(length(n),1);
for i = 1 : length(n)
    % Face's tangential vector
    t(i,:) = cross( cross(k,n(i,:)) , n(i,:) ) / norm( cross( k,n(i,:) ) );
    % Impact Angle, the angle between relative velocity vector and face's
    % tangential vector
    theta_impact(i) = acos(  dot( v_r_norm,t(i,:)' ) / norm(v_r_norm) / norm(t(i,:)')  ) ;
end


% Calculate the thetaOUT based on S.D. Raducan et.al work: 
% Ejecta distribution and momentum transfer from oblique impacts on
% asteroid surfaces
gamma = zeros(length(n),1);
for i = 1 : length(theta_impact)
    if theta_impact(i) > pi/2
        x = theta_impact(i) - pi/2;
    else
        x = theta_impact(i);
    end
    x = x*180/pi;
    gamma(i) =1 +  (beta-1)*tan(x*pi/180)*tan((-0.0052*x^2+1.3629*x-80.871)*pi/180);
end


% Ready to calculate delta_v to PHA
norm_delta_v_A = zeros(length(gamma),1);
delta_v = zeros(length(gamma),3);
% delta_v_mass_point = zeros(length(beta),3);
for i = 1 : length(gamma)
    % delta_v
    delta_v(i,1:3) = m/M * (...
        beta*dot( v_r,n(i,:)' )*n(i,:)' + ...
        gamma(i)*dot(v_r,t(i,:)')*t(i,:)' ...
        );
    norm_delta_v_A(i,1) = abs( dot(delta_v(i,1:3),v_PHA) / norm(v_PHA) );
end
% Find the index that generated largest proj_delta_v to v_PHA
[~,b] = max( norm_delta_v_A);
% This delta_v is the best impact point (BIP) described in the paper
best_delta_v = delta_v(b,1:3);
% Store the attitude
angles_output(1,:) = angles(:,1);
% -------------------------------------------------------------------------





% -------------------- COG strategy (Center of Mass) ----------------------
% Calculate the faces that irradiated_vector pass though its center of mass
center_of_mass = mean(vertices);
ray_origin = center_of_mass;
ray_end = ray_origin - irradiated_vector' * 10;
effective_faces = faces(irradiated_face_indices, :);
vert1 = rotated_vertices(effective_faces(:,1),:);
vert2 = rotated_vertices(effective_faces(:,2),:);
vert3 = rotated_vertices(effective_faces(:,3),:);
% Check if the vector passes through the center of mass of the face
intersect = TriangleRayIntersection(ray_origin, ray_end, vert1, vert2, vert3);
face_through_COM_indices = intersect == 1;

% Calculate the tangential_vector_t & theta_impact of COM
n_COM = -face_normals(irradiated_face_indices,:);  %normal components of face
n_COM = n_COM(face_through_COM_indices,:);
k = [0,0,1];
try
    t_COM = cross(cross(k, n_COM), n_COM) / norm(cross(k, n_COM));
catch ME
    disp('Error in cross product calculation:');
    disp(ME.message);
    disp('Size of k:');
    disp(size(k));
    disp('Size of n_COM:');
    disp(size(n_COM));
    disp('Values of k and n_COM:');
    disp('k =');
    disp(k);
    disp('n_COM =');
    disp(n_COM);
    disp('This error is caused by during the TriangleRayIntersection (in get_delta_v_from_momentum.m, line: 140) to find COG impact surfaces.')
    disp('The 3D model of PHA is not fine enough, the ray failed to intersect a surface (it went through edege)')
    disp('This is a small propability case, rerun the code should be fine. Or can switch to a more detailed 3D PHA model')
    rethrow(ME); % 或者使用 error 自定义错误
end
t_COM = cross( cross(k,n_COM) , n_COM ) / norm( cross( k,n_COM ) );
theta_impact_COM = acos(  dot( v_r_norm,t_COM' ) / norm(v_r_norm) / norm(t_COM')  ) ;


    if theta_impact_COM > pi/2
        x_COG = theta_impact_COM - pi/2;
    else
        x_COG = theta_impact_COM;
    end
    x_COG = x_COG*180/pi;
    gamma_COM =1 +  (beta-1)*tan(x_COG*pi/180)*tan((-0.0052*x_COG^2+1.3629*x_COG-80.871)*pi/180);

% delta_v generated through COG impact
delta_v_COM = m/M * (...
        beta*dot( v_r,n_COM' )*n_COM' + ...
        gamma_COM*dot(v_r,t_COM')*t_COM' ...
        );
% -------------------------------------------------------------------------

end
