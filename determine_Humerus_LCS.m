function [ T, R, ori, success_status ] = determine_Humerus_LCS( marker_coords )

% Assuming that motion capture marker data for 4 markers (!!)
% is given in a variable matrix marker_coords (4x3, 
% i.e. n_markers x n_coordinates --> only data for 1 timepoint)
% the function calculates a 4x4 transformation matrix 
% that maps from the motion capture to the local probe coordinate system.
%
% The assumed design of the marker cluster data is such that there is
% 1st row: marker Lax1 (Long axis, point 1) 
% 2nd row: marker Lax2 (Long axis, point 2) 
% 3rd row: marker Sax1 (Short axis, point 1) 
% (coordinates of the 5th marker point are not used here)

% initialise return values
success_status = 0;
T   = [];
R   = [];
ori = [];

% assign specific variables to each of the points acc. to the above defs
Lax1 = marker_coords(2,:);
Lax2 = marker_coords(3,:);

Sax1 = marker_coords(1,:);

% define two unit direction vectors to establish a reference plane
u = Lax1-Lax2;
u = u/norm(u);

ori = Lax1 + dot(Sax1-Lax1,Lax2-Lax1) / dot(Lax2-Lax1,Lax2-Lax1) * (Lax2-Lax1);

v = Sax1-ori;
v = v/norm(v);

% calculate the normal vector to the reference plane
n = cross(u,v);
n = n/norm(n);

% assign/re-name initial unit direction vectors to establish R matrix
ex = u;
ey = v;
ez = n;

% to establish 3 mutually perpendicular axes required for R, 
% we can use n (ez), and also u (ex, from the longer axis) "as is"; 
% however, initially v (ey) is not exactly orthogonal to u (ex),
% so we ensure it is now:
ey = cross(ez, ex);
ey = ey/norm(ey);

% we have the mutually perpendicular unit direction vectors to serve as
% basis vectors to define a 3x3 rotation matrix now!
R = [ ex; ey; ez ];

% create a 4x4 transformation matrix that maps to the local CS 
T = create_4by4_tmatrix( R, ori );
success_status = 1;


return
end