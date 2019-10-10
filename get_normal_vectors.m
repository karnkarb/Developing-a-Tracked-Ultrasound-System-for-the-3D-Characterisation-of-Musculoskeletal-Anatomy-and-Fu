function [ u, v ] = get_normal_vectors( w )

% given a 3D vector w, find two vectors normal to w such that:
% u x v = w
% v x w = u
% w x u = v

% make sure w is indeed a unit vector
w = w/norm(w);

% find a vector that is normal to w
[ ~, min_idx ] = min(w);
[ ~, max_idx ] = max(w);
non_minmax_idx = setdiff([1,2,3], [min_idx, max_idx]);

% check for special case: min = non_min = 0 !!!
eps_zero = 1e-9;

if ( abs(min(w(min_idx)))        < eps_zero && ...
     abs(min(w(non_minmax_idx))) < eps_zero     )

    u_tmp(min_idx)        = 0.0;
    u_tmp(non_minmax_idx) = 1.0;
    u_tmp(max_idx)        = 0.0;     
    u_tmp = u_tmp/norm(u_tmp);

    v = cross(w,u_tmp);
    v = v/norm(v);
    
    u = cross(v,w);
    u = u/norm(u);
else
    u(min_idx) = 0.0;
    u(non_minmax_idx) = -w(max_idx);
    u(max_idx)        =  w(non_minmax_idx);
    u = u/norm(u);

    v = cross(w,u);
    v = v/norm(v);
end

% test: e.g.: R [ u; v; w ]; det(r) == 1 ???

return

end