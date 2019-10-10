function [ OPA_Marker_Position, Rotation, Translation ] ...
   = OPA ( Source_Marker_Position, Target_Marker_Position )

%-----------------------------------------------------------------------
%
% This program computes the Ordinary Procrustes Analysis (OPA)
% for marker data.
%
% Reference : I. Söderkvist, P.A. Wedin, Determining the movements
%             of the skeleton using well-configured markers.
%             J. Biomech. 26, 1473-7 (1993)
%
%-----------------------------------------------------------------------
%   v2.1 - 08-02-2015
%
%   Written by:
%
%   Rainald Ehrig
%   Zuse Institute Berlin (ZIB)
%   Takustrasse 7
%   D-14195 Berlin-Dahlem
%   Phone : +49-30-84185-0
%   Fax   : +49-30-84185-125
%   E-mail: ehrig@zib.de
%   URL   : http://www.zib.de
%
%   change log:
%   v2.1 revised centring of the markers 
%        08-02-2015 Markus Heller m.o.heller@soton.ac.uk 
%   v2.2 vectorised loops 
%        12-03-2015 Markus Heller m.o.heller@soton.ac.uk 
%
%-----------------------------------------------------------------------
%
% Implemented without slow features like repmat or mean!
%
%-----------------------------------------------------------------------
%
% Centre the configurations
%
%-----------------------------------------------------------------------

[ ~, Number_of_Markers ] = size ( Source_Marker_Position );

% cS = zeros ( 3, 1 );
% cT = zeros ( 3, 1 );

cS = sum(Source_Marker_Position, 2) / Number_of_Markers;
cT = sum(Target_Marker_Position, 2) / Number_of_Markers;

% centS = zeros ( 3, Number_of_Markers );
% centT = zeros ( 3, Number_of_Markers );

centS = bsxfun(@minus, Source_Marker_Position, cS);
centT = bsxfun(@minus, Target_Marker_Position, cT);

% for i = 1 : Number_of_Markers;
% 
%     centS(:,i) = Source_Marker_Position(:,i) - cS;
%     centT(:,i) = Target_Marker_Position(:,i) - cT;
%     
% end;

%-----------------------------------------------------------------------
%
% Compute SVD
%
%-----------------------------------------------------------------------

[ U, ~, V ] = svd ( centT * centS' );

%----------------------------------------------------------------------
%
% Compute uv' and check determinant from uv', if negative recompute uv'
%
%-----------------------------------------------------------------------

Rotation = U * V';

if ( det  ( Rotation )  < 0 )
    
   V(:,3) = - V(:,3);
   Rotation = U * V';

end

%-----------------------------------------------------------------------
%
% Determine translation vector
%
%-----------------------------------------------------------------------

Translation = cT - Rotation * cS;

%-----------------------------------------------------------------------
%
% Compute transformed source configuration
%
%-----------------------------------------------------------------------
                  
OPA_Marker_Position = Rotation * Source_Marker_Position;

OPA_Marker_Position = bsxfun(@plus, OPA_Marker_Position, Translation);

% for i = 1 : Number_of_Markers;
%     OPA_Marker_Position(:,i) = OPA_Marker_Position(:,i) + Translation;
% end;

%-----------------------------------------------------------------------

return;
end