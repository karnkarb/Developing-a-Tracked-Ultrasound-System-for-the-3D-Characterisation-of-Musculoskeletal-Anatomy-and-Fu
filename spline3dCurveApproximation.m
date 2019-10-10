function [ spline_pts, spline_derivs, FS_frame, spline_result ] = spline3dCurveApproximation ( pts, n_knots, n_splev, curve_type ) 

% approximates a 3d curve (open or closed) using splines
% the larger n_knots is, the closer will the spline follow the data points

% uses SPLINEFIT by Jonas Lundgren
% https://uk.mathworks.com/matlabcentral/fileexchange/13812-splinefit
% Set number of spline evaluations - define as parameter?

%--------------------------------------------------------------------------
% inputs:
% =======
% pts         nx3,            3d point coordinates to be fitted
% n_knots     scalar value,   number of knots used for the spline - the more
%                             knots, the closer the spline will follow the
%                             original data (i.e. less smoothing)
% n_splev     scalar value,   number of spline evaluations
%                             how densely do we want the spline to be evaluated?
%                             larger numbers will mean a more accurate visual 
%                             representation (order: a few hundreds to 1000 - 4000)
% curve_typ   scalar value,   0: open; 1: closed curve
%
% output:
% =======
% spline_pts    mx3,            3d coordinates of the spline
% spline_deriv  mx3             1st derivative of spline coords
% spline_result structure       actual spline data
%               .xspl, .yspl, z.spl
%               .xspl_deriv, .yspl_deriv, .zspl_deriv
%--------------------------------------------------------------------------
% v1.0 28.02.2018 (c) M.O. Heller markus.o.heller@soton.ac.uk
% v1.1 27.08.2018     added 1st derivative, return actual spline info
%--------------------------------------------------------------------------

n_pts = size(pts, 1);

%--------------------------------------------------------------------------
% proceed with the actual spline fit
%--------------------------------------------------------------------------

% sort points into vectors for x, y and z coordinates
x(1:n_pts) = pts(:, 1); 
y(1:n_pts) = pts(:, 2); 
z(1:n_pts) = pts(:, 3);

% calculate an appropriate parametrisation based on chord length
t = cumsum ( sqrt([0,diff(x)].^2 + [0,diff(y)].^2 + [0,diff(z)].^2 ) );

% perform a spline fit considering that we either have a closed or
% an open curve to fit here
if ( curve_type == 0 )  % open curve
   splinex = splinefit ( t, x, n_knots );
   spliney = splinefit ( t, y, n_knots );
   splinez = splinefit ( t, z, n_knots );
else                    % closed curve
   splinex = splinefit ( t, x, n_knots, 'p' );
   spliney = splinefit ( t, y, n_knots, 'p' );
   splinez = splinefit ( t, z, n_knots, 'p' );
end

tt = linspace ( t(1), t(end), n_splev );
        
% at these points we evaluate the splines
xspl = ppval ( splinex, tt );
yspl = ppval ( spliney, tt );
zspl = ppval ( splinez, tt );

% prepare data to return
spline_pts   = ([ xspl; yspl; zspl ])';


% let's also calculate the 1st derivative
splinex_deriv1 = ppdiff ( splinex, 1 );
spliney_deriv1 = ppdiff ( spliney, 1 );
splinez_deriv1 = ppdiff ( splinez, 1 );

splinex_deriv2 = ppdiff ( splinex, 2 );
spliney_deriv2 = ppdiff ( spliney, 2 );
splinez_deriv2 = ppdiff ( splinez, 2 );

xspl_deriv1 = ppval ( splinex_deriv1, tt );
yspl_deriv1 = ppval ( spliney_deriv1, tt );
zspl_deriv1 = ppval ( splinez_deriv1, tt );

xspl_deriv2 = ppval ( splinex_deriv2, tt );
yspl_deriv2 = ppval ( spliney_deriv2, tt );
zspl_deriv2 = ppval ( splinez_deriv2, tt );

spline_deriv1 = ([ xspl_deriv1; yspl_deriv1; zspl_deriv1 ])';
spline_deriv2 = ([ xspl_deriv2; yspl_deriv2; zspl_deriv2 ])';

% let's try to calculate a local Frenet-Serret frame / CS -----------------
% following Farouki 2016: Rational rotation-minimizing frames - Recent
% advances and open problems.
%--------------------------------------------------------------------------
T = spline_deriv1./vecnorm(spline_deriv1, 2, 2);  % tangent T
tmp_b = cross(spline_deriv1, spline_deriv2);
B = tmp_b./vecnorm(tmp_b, 2, 2); % binormal B
P = cross(B, T);
P = P./vecnorm(P, 2, 2); % principal normal P

FS_frame.T = T;
FS_frame.P = P;
FS_frame.B = B;


% also have a go at calculating curvature and torsion
tmp_kap1 = vecnorm(cross(spline_deriv1, spline_deriv2), 2, 2);
tmp_kap2 = vecnorm(spline_deriv1, 2, 2).^3;
kappa = tmp_kap1./tmp_kap2;

FS_frame.kappa = kappa;
FS_frame.R     = 1./kappa;


% return the data set
spline_derivs.first  = spline_deriv1;
spline_derivs.second = spline_deriv2;


spline_result.xspl = xspl;
spline_result.yspl = yspl;
spline_result.zspl = zspl;

spline_result.xspl_deriv1 = splinex_deriv1;
spline_result.yspl_deriv1 = spliney_deriv1;
spline_result.zspl_deriv1 = splinez_deriv1;

spline_result.xspl_deriv2 = splinex_deriv2;
spline_result.yspl_deriv2 = spliney_deriv2;
spline_result.zspl_deriv2 = splinez_deriv2;

return
end