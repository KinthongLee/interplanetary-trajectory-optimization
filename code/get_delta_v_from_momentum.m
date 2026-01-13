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
t = zeros(length(n),3);  % tangential components of face
theta_impact = zeros(length(n),1);
for i = 1:size(n,1)

    vr = v_r_norm(:).';     % ensure row vector
    ni = n(i,:);

    cn = cross(ni,vr);

    % For some special case when vr nearly parallel wiht ni, there are
    % infinity of possible t. Simply give a k2 vector
    if norm(cn) < 1e-10
        k2 = [0,0,1];
        if abs(dot(k2,ni)) > 0.99
            k2 = [1,0,0];
        end
        cn = cross(k2, ni);
    end
    
    t(i,:) = cross(cn, ni) / norm(cn);

    % Make sure t & vr are same direction
    if dot(t(i,:), vr) < 0
        t(i,:) = -t(i,:);
    end

    % Impact Angle
    theta_impact(i) = acos( dot(vr, t(i,:)) ...
                           /(norm(vr)*norm(t(i,:))) );
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
center_of_mass = mean(rotated_vertices, 1);
ray_origin = center_of_mass;
ray_end = ray_origin - irradiated_vector' * 10;

effective_faces = faces(irradiated_face_indices, :);
vert1 = rotated_vertices(effective_faces(:,1),:);
vert2 = rotated_vertices(effective_faces(:,2),:);
vert3 = rotated_vertices(effective_faces(:,3),:);

% Check if the vector passes through the center of mass of the face
intersect = TriangleRayIntersection(ray_origin, ray_end, vert1, vert2, vert3);
face_through_COM_indices = intersect == 1;

% Face normal of COM impact surface (should be exactly one face in normal cases)
n_COM_all = -face_normals(irradiated_face_indices,:);   % outward normal facing spacecraft
n_COM = n_COM_all(face_through_COM_indices,:);

% Robust guard: ray intersection should hit exactly one triangle face.
% If multiple faces are hit, pick the one most "facing" the incoming direction.
if isempty(n_COM)
    error(['COG surface not found: TriangleRayIntersection returned no hit. ', ...
           'Likely the ray went through an edge; use a finer mesh or rerun it.',...
           'This is a small probability case, because the default Apophis light-curve model is not fine enough',...
           'just rerun the code could fix it']);
end
if size(n_COM,1) > 1
    % choose face with most negative dot(face_normal, irradiated_vector) i.e. most facing spacecraft
    dp = dot(n_COM, repmat(irradiated_vector(:).', size(n_COM,1), 1), 2);
    [~, idxBest] = min(dp);     % more negative -> better facing
    n_COM = n_COM(idxBest,:);
end

% Ensure n_COM is unit row vector
n_COM = n_COM(:).';
n_COM = n_COM / norm(n_COM);

% --------- Construct tangential vector t_COM to be coplanar with v_r and n_COM ----------
% Define t_COM as the projection of v_r onto the tangent plane of the face:
% t_COM ∝ v_r_norm - (v_r_norm·n_COM) n_COM
% This guarantees: t_COM ⟂ n_COM, and {n_COM, t_COM, v_r_norm} are coplanar.
vr = v_r_norm(:).';   % row, unit

cn = cross(n_COM, vr);          % this is perpendicular to both v_r and n
if norm(cn) < 1e-12
    % Degenerate: v_r parallel to n (near-normal impact), tangential direction is not unique.
    % Pick a deterministic fallback axis not parallel to n_COM to define a stable tangent.
    k2 = [0,0,1];
    if abs(dot(k2, n_COM)) > 0.99
        k2 = [1,0,0];
    end
    cn = cross(k2, n_COM);
end

t_COM = cross(cn, n_COM);
t_COM = t_COM / norm(t_COM);   % unit row vector

% % Optional: enforce a consistent sign so that t_COM points "with" the tangential component of v_r
if dot(t_COM, vr) < 0
    t_COM = -t_COM;
end

% --------- Impact angle definition ----------
theta_impact_COM = acos( dot(vr, t_COM) /norm(vr) / norm(t_COM) );


% Your existing gamma mapping uses x in degrees with a piecewise for theta > pi/2
if theta_impact_COM > pi/2
    x_COG = theta_impact_COM - pi/2;
else
    x_COG = theta_impact_COM;
end
x_COG = x_COG * 180/pi;

gamma_COM = 1 + (beta-1) * tan(x_COG*pi/180) * tan( (-0.0052*x_COG^2 + 1.3629*x_COG - 80.871)*pi/180 );

% delta_v generated through COG impact (row vector 1x3)
delta_v_COM = (m/M) * ( ...
    beta    * dot(v_r, n_COM') * n_COM' + ...
    gamma_COM * dot(v_r, t_COM') * t_COM' ...
    );

% -------------------------------------------------------------------------


end
