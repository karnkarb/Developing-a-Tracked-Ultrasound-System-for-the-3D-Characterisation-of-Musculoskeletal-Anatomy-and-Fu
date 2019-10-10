function [ GPA_Marker_Position, conv ] = GPA ( Marker_Position )

%-----------------------------------------------------------------------
%
% This program computes the Generalized Procrustes Analysis (GPA)  
% without scaling, i.e. the Optimal Common Shape for a set of markers.
%
% References: I.L. Dryden, K. Mardia, Statistical Shape Analysis,
%             Wiley, 2002
%
%-----------------------------------------------------------------------
%   v2.0 21-03-2014
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
%-----------------------------------------------------------------------
% 
% GPA requires as input a set of marker positions in the 3-dimensional
% field Marker_Position where the size of the array dimensions are:
% 3 (the x, y, z coordinates), number of markers, number of time frames.
%
% The output of GPA is a set of corrected marker positions of the same
% size as the input data set. These marker positions have a fixed 
% geometry, i.e. all marker-marker distances the same for each time
% frame. 
%
%-----------------------------------------------------------------------
%
% Project all frames onto first frame
%
%-----------------------------------------------------------------------

[ ~, Number_of_Markers, Number_of_Frames ] = size ( Marker_Position );

GPA_Marker_Position = zeros ( 3, Number_of_Markers, Number_of_Frames );

for i = 1 : Number_of_Frames
   [ GPA_Marker_Position(:,:,i), ~, ~ ] ...
      = OPA ( Marker_Position(:,:,i), Marker_Position(:,:,1) );
end

%-----------------------------------------------------------------------
%
% Centre all configurations
%
%-----------------------------------------------------------------------

Mean_Marker_Position = mean ( GPA_Marker_Position, 2 );

for j = 1 : Number_of_Markers
   GPA_Marker_Position(:,j,:) ...
      = GPA_Marker_Position (:,j,:) - Mean_Marker_Position(:,1,:);
end

%-----------------------------------------------------------------------
% 
% Set initial value for the residual
%
%-----------------------------------------------------------------------

res1 = 1.0d20;
conv = 0;

%-----------------------------------------------------------------------
% 
% Compute mean configuration
%
%-----------------------------------------------------------------------

Mean_Configuration = mean ( GPA_Marker_Position, 3 );

%-----------------------------------------------------------------------
%
% Begin of main GPA iteration
%
%-----------------------------------------------------------------------

for iGPA = 1 : 10
    
%-----------------------------------------------------------------------
%
% Compute rotations
%
%-----------------------------------------------------------------------

   for i = 1 : Number_of_Frames
      [ GPA_Marker_Position(:,:,i), ~, ~ ] ...
         = OPA ( GPA_Marker_Position(:,:,i), Mean_Configuration );
   end
   
   Mean_Configuration = mean ( GPA_Marker_Position, 3 );

   res2 = 0;
   
   for k = 1 : Number_of_Frames
      res2 = res2 + norm ( GPA_Marker_Position(:,:,k) - Mean_Configuration, 'fro' ) ^ 2;
   end

   if ( res2 / res1 > 0.9999 ) 
     
      conv = 1;
      break;

   end
      
   res1 = res2;

end

%-----------------------------------------------------------------------
%
% Print warning if no convergence achieved in 10 iterations
%
%-----------------------------------------------------------------------

if conv == 0
   fprintf ( '\n More than 10 iterations needed for convergence of GPA, check data!\n' );
end

%-----------------------------------------------------------------------
%
% Back-projection of the mean configuration to the original positions
%
%-----------------------------------------------------------------------

for i = 1 : Number_of_Frames
   [ GPA_Marker_Position(:,:,i), ~, ~ ] ...
      = OPA ( Mean_Configuration, Marker_Position(:,:,i) );
end

%-----------------------------------------------------------------------

return
end 