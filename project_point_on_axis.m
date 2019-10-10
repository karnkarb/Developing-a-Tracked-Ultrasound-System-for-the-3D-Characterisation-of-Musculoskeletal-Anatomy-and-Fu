function [ q, t ] = project_point_on_axis( p, r, v )
%project point on axis
% project point p on the axis given by point r and direction vector v, return result as q
% MHeller 12.01.2014

% the axis is parameterized as: r + t  * v

% determine the vector from point r (on the axis) to point p, which needs to be projected onto the axis
rp = p-r;

% determine the projection of point p onto the axis - the point is found at:
% t = [(rp) . (v)] / |v|^2
v_squared = dot(v,v);    
t = dot(rp,v)/v_squared;  
q = r + t * v;
  
end