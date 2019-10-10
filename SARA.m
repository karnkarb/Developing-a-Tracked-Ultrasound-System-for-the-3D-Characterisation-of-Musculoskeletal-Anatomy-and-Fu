%function [ Rotation_1, Rotation_2, ...
%           Translation_1, Translation_2, ...
%           SVD_Matrix, RHS, ...
%           U, V, ...
function [ Singular_Values, ...% Reduced_Singular_Values, ...
           Residual, Scaled_Residual, ...
           Angular_Variance, Range_Of_Motion, ...
           Local_Centres, Local_Axis_Points, Local_Axis_Directions, ...
           Global_Centres_1, Global_Centres_2, Mean_Global_Centres, ...
           Global_Axis_Points_1, Global_Axis_Directions_1, ... 
           Global_Axis_Points_2, Global_Axis_Directions_2, ...
           Global_Mean_Axis_Points, Global_Mean_Axis_Directions ] ...
    = SARA ( Marker_Position_1, Marker_Position_2 )

%-----------------------------------------------------------------------
%
% This program computes joint centre and axis positions
%
%-----------------------------------------------------------------------
%
%     Written by:
%     R. Ehrig
%     Zuse Institute Berlin (ZIB)
%     Takustrasse 7
%     D-14195 Berlin-Dahlem
%     phone : +49-30-84185-0
%     e-mail: ehrig@zib.de
%     URL   : http://www.zib.de
%
%-----------------------------------------------------------------------
% 
% Input:
% ======
%
% Marker_Position_1 : marker positions of segment 1 
% in a 3 dimensional matrix with dimensions:
% 3 (the x, y, z coordinates), number of markers, number of time frames
% (as in GPA.m)
%
% Marker_Position_2 : marker positions of segment 2
% in a 3 dimension matrix as above
%
% Output:
% =======
%
% Singular_Values : the singular values of the SCoRE matrix, vector of 
%     length 6
%
% Residual : the residual of the overdetermined linear SCoRE system
%
% Scaled_Residual : the scaled residual (residual per time frame)
%
% Angular_Variance : the angular variance of the rotations around the
%     principal rotation axes, vector of length 3
%
% Range_Of_Motion : the range of motion around the three principal axes,
%     vector of length 3
%
% Local_Centres : the centre position of the local representations of the 
%    joint centre, vector of length 6, Local_Centres(1:3) is the
%    representation of the local centre in segment 1, Local_Centres(4:6)
%    the representation of the local centre in segment 2
%
% Local_Axis_Points : the points on the three principal rotation joint
%    axes nearest to the local origin, matrix of  dimensions 6 x 3,
%    Local_Axis_Points(1:3,i) is the axis point on the i-th principal
%    rotation axis in local coordinates of segment 1, 
%    Local_Axis_Points(4:6,i) correspondingly those in segment 2
%
% Local_Axis_Directions : the normalized direction vectors of the three
%    principal rotation joint axes, matrix of dimensions 6 x 3, 
%    Local_Axis_Direction(1:3,i) is the direction vector of the i-th 
%    principal rotation axis in local coordinates of segment 1,
%    Local_Axis_Direction(4:6,i) correspongly those of segment 2
%
% Global_Centres_1 : the global time dependent centre representation of 
%    the joint centre associated to segment 1, matrix of dimension
%    3 x Number_of_Frames, thus Global_Centres(1:3,i) describes the
%    joint centre position w.r.t. segment 1 at time frame i
%
% Global_Centres_2 : the global time dependent centre representation of 
%    the joint centre associated to segment 2, matrix of dimension
%    3 x Number_of_Frames, as above for segment 2
%
% Mean_Global_Centres : the mean between the two time dependent 
%    representations of the global joint centre, matrix of dimension
%    3 x Number_of_Frames, provably the statistically optimal
%    estimate for the joint centre postion in global coordinates
%
% Global_Axis_Points_1 : the time dependent points on the three principal 
%    rotation joint axes associated to segment 1 nearest to the global 
%    origin, matrix of  dimensions 3 x 3 x Number_of_Frames, thus 
%    Global_Axis_Points_1(1:3,i,j) is the axis point on the i-th principal
%    rotation axis w.r.t. segment 1 at time frame j in global coordinates
%
% Global_Axis_Directions_1 : the time dependent normalized direction 
%    vectors of the three principal rotation joint axes associated to 
%    segment 1, matrix of  dimensions 3 x 3 x Number_of_Frames, thus 
%    Global_Axis_Directions_1(1:3,i,j) is the direction vector of the i-th 
%    principal rotation axis w.r.t. segment 1 at time frame j in global 
%    coordinates,
% 
% Global_Axis_Points_2 : the time dependent points on the three principal 
%    rotation joint axes associated to segment 2 nearest to the global 
%    origin, matrix of dimensions 3 x 3 x Number_of_Frames, as above for 
%    segment 1
%
% Global_Axis_Directions_2 : the time dependent normalized direction 
%    vectors of the three principal rotation joint axes associated to
%    segment 2, matrix of dimensions 3 x 3 x Number_of_Frames, as above 
%    for segment 1
%
% Mean_Global_Axis_Points : the time dependent points on the mean first
%    principal rotation joint axis nearest to the global origin, matrix 
%    of  dimensions 3 x Number_of_Frames
%
% Mean_Global_Axis_Directions : the time dependent normalized direction 
%    vectors of the first principal rotation joint axis, matrix of 
%    dimensions 3 x Number_of_Frames
%
%-----------------------------------------------------------------------


% Determine the number of frames
[ ~, ~, Number_of_Frames ] = size ( Marker_Position_1 );

%-----------------------------------------------------------------------
% Construct local coordinate systems in both segments,
% by determining rotations and translations
%-----------------------------------------------------------------------
Rotation_1 = zeros ( 3, 3, Number_of_Frames ); % segment 1
Rotation_2 = zeros ( 3, 3, Number_of_Frames ); % segment 2

% First axis direction vectors
Rotation_1(:,1,:) = Marker_Position_1(:,2,:) - Marker_Position_1(:,1,:);
Rotation_2(:,1,:) = Marker_Position_2(:,2,:) - Marker_Position_2(:,1,:);

% Second axis direction vectors
Rotation_1(:,2,:) = Marker_Position_1(:,3,:) - Marker_Position_1(:,1,:);
Rotation_2(:,2,:) = Marker_Position_2(:,3,:) - Marker_Position_2(:,1,:);

% Normalize, orthogonalize, and compute cross products (third axis direction vectors)
for i = 1 : Number_of_Frames
    % normalise the direction vectors of the 1st axis for segment 1 and 2
    Rotation_1(:,1,i) = Rotation_1(:,1,i) / norm ( Rotation_1(:,1,i) );
    Rotation_2(:,1,i) = Rotation_2(:,1,i) / norm ( Rotation_2(:,1,i) );
    
    % 
    Rotation_1(:,2,i) = Rotation_1(:,2,i) ... 
       - ( Rotation_1(:,1,i) * Rotation_1(:,2,i)' ) * Rotation_1(:,1,i);
    Rotation_2(:,2,i) = Rotation_2(:,2,i) ... 
       - ( Rotation_2(:,1,i) * Rotation_2(:,2,i)' ) * Rotation_2(:,1,i);
    
    % normalise the direction vectors of the 2nd axis for segment 1 and 2 
    Rotation_1(:,2,i) = Rotation_1(:,2,i) / norm ( Rotation_1(:,2,i) );
    Rotation_2(:,2,i) = Rotation_2(:,2,i) / norm ( Rotation_2(:,2,i) );
    
    % determine the direction vectors of the 3rd axis using the cross product
    Rotation_1(:,3,i) = cross ( Rotation_1(:,1,i),  Rotation_1(:,2,i) );
    Rotation_2(:,3,i) = cross ( Rotation_2(:,1,i),  Rotation_2(:,2,i) );
end

% Detemine the translation as 
Translation_1 = squeeze ( mean ( Marker_Position_1(:,1:3,:), 2 ) );
Translation_2 = squeeze ( mean ( Marker_Position_2(:,1:3,:), 2 ) );


% disp('dimensions of Rotation_1/2 are: ')
% size(Rotation_1)
% disp('dimensions of Translation_1/2 are: ')
% size(Translation_1)


% Construct local coordinate systems ---------------------------------- END


% Set-up the matrix for calculating the SVD -------------------------------
SVD_Matrix = zeros ( 3 * Number_of_Frames, 6 );

for i = 1 : Number_of_Frames
    SVD_Matrix(3*i-2:3*i,:) = ...
        horzcat ( Rotation_1(:,:,i), -Rotation_2(:,:,i) );
end

% disp('dimensions of SVD_Matrix are: ')
% size(SVD_Matrix)
% Set-up the matrix for calculating the SVD --------------------------- END


% Define the Right Hand Side (RHS)
RHS = reshape ( Translation_2 - Translation_1, 3 * Number_of_Frames, 1 );

% disp('dimensions of RHS are: ')
% size(RHS)


% Perform the SVD and assign singular values ------------------------------
[ U, Singular_Values, V ] = svd ( SVD_Matrix, 0  ); 

Singular_Values = diag ( Singular_Values );
% Perform the SVD and assign singular values -------------------------- END


% disp('dimensions of U are: ')
% size(U)
% disp('dimensions of Singular_Values are: ')
% size(Singular_Values)
% disp('dimensions of V are: ')
% size(V)



% Determine the reduced singular values
Reduced_Singular_Values = ...
    Singular_Values(1:3).^2 / Number_of_Frames - ones ( 3, 1 );

% disp('dimensions of Reduced_Singular_Values are: ')
% size(Reduced_Singular_Values)


% Determine angular variances         
Angular_Variance = zeros ( 3, 1 );

Angular_Variance(1) = 1.0 - sqrt (   Reduced_Singular_Values(2) ...
                                   * Reduced_Singular_Values(3) ...
                                   / Reduced_Singular_Values(1) );
Angular_Variance(2) = 1.0 - sqrt (   Reduced_Singular_Values(1) ...
                                   * Reduced_Singular_Values(3) ...
                                   / Reduced_Singular_Values(2) );
Angular_Variance(3) = 1.0 - sqrt (   Reduced_Singular_Values(1) ...
                                   * Reduced_Singular_Values(2) ...
                                   / Reduced_Singular_Values(3) );

% Determine ranges of motion (ROM)
Range_Of_Motion = RangeOfMotion ( Angular_Variance );


% Compute the norm minimal solution
Local_Centres = zeros ( 6, 1 );

% for i = 1 : 6
%     Local_Centres(:) = ...
%         Local_Centres(:) + V(:,i) * ( RHS' * U(:,i) ) / Singular_Values(i);
% end

for i = 1 : 6
    Local_Centres = Local_Centres + V(:,i) * ( RHS' * U(:,i) ) / Singular_Values(i);
end

% Compute the residual
Residual        = norm ( SVD_Matrix * Local_Centres - RHS );
Scaled_Residual = Residual / sqrt ( Number_of_Frames );


% Compute the local axes from the right singular vectors
Local_Axis_Directions = zeros ( 6, 3 );
Local_Axis_Points     = zeros ( 6, 3 );

for i = 1 : 3    
    Local_Axis_Directions(1:3,i) =   V(1:3,i) / norm ( V(1:3,i) );
    Local_Axis_Directions(4:6,i) = - V(4:6,i) / norm ( V(4:6,i) );
end

% Contruct axis points nearest to local origin

for i = 1 : 3    
    Local_Axis_Points(1:3,i) = Local_Centres(1:3) ...
        -   ( Local_Centres(1:3)' * Local_Axis_Directions(1:3,i) ) ...
          * Local_Axis_Directions(1:3,i);
    Local_Axis_Points(4:6,i) = Local_Centres(4:6) ...
        -   ( Local_Centres(4:6)' * Local_Axis_Directions(4:6,i) ) ...
          * Local_Axis_Directions(4:6,i);
end

% Backtransformation to global coordinates --------------------------------

% Centres
Global_Centres_1 = zeros ( 3, Number_of_Frames );
Global_Centres_2 = zeros ( 3, Number_of_Frames );

for i = 1 : Number_of_Frames    
    Global_Centres_1(:,i) =   Translation_1(:,i) ...
                            + Rotation_1(:,:,i) * Local_Centres(1:3);
    Global_Centres_2(:,i) =   Translation_2(:,i) ...
                            + Rotation_2(:,:,i) * Local_Centres(4:6);    
end

Mean_Global_Centres = 0.5 * ( Global_Centres_1 + Global_Centres_2 );

% Axes
Global_Axis_Points_1     = zeros ( 3, 3, Number_of_Frames );
Global_Axis_Points_2     = zeros ( 3, 3, Number_of_Frames );
Global_Axis_Directions_1 = zeros ( 3, 3, Number_of_Frames );
Global_Axis_Directions_2 = zeros ( 3, 3, Number_of_Frames );

for i = 1 : 3    
    for j = 1 : Number_of_Frames        
        Global_Axis_Points_1(:,i,j) = Translation_1(:,j) ...
            + Rotation_1(:,:,j) * Local_Axis_Points(1:3,i);
        Global_Axis_Points_2(:,i,j) = Translation_2(:,j) ... 
            + Rotation_2(:,:,j) * Local_Axis_Points(4:6,i);
    
        Global_Axis_Directions_1(:,i,j) = ...
            Rotation_1(:,:,j) * Local_Axis_Directions(1:3,i);
        Global_Axis_Directions_2(:,i,j) = ...
            Rotation_2(:,:,j) * Local_Axis_Directions(4:6,i);
    end
end


% Construct axis points nearest to global origin 
for i = 1 : 3    
   for j = 1 : Number_of_Frames        
      Global_Axis_Points_1(:,i,j) = Global_Axis_Points_1(:,i,j) ...
 -   ( Global_Axis_Points_1(:,i,j)' * Global_Axis_Directions_1(:,i,j) ) ...
   * Global_Axis_Directions_1(:,i,j);
      Global_Axis_Points_2(:,i,j) = Global_Axis_Points_2(:,i,j) ...
 -   ( Global_Axis_Points_2(:,i,j)' * Global_Axis_Directions_2(:,i,j) ) ...
   * Global_Axis_Directions_2(:,i,j);        
   end
end

% Compute mean global axes

Global_Mean_Axis_Points     = zeros ( 3, Number_of_Frames );
Global_Mean_Axis_Directions = zeros ( 3, Number_of_Frames );

for i = 1 : Number_of_Frames
   [ Global_Mean_Axis_Points(:,i), Global_Mean_Axis_Directions(:,i) ] ...
      = MeanAxis ( Global_Axis_Points_1(:,1,i), Global_Axis_Directions_1(:,1,i), ...
                   Global_Axis_Points_2(:,1,i), Global_Axis_Directions_2(:,1,i) );
               
%                if ( i == 1 )
%                    disp('dimensions of Global_Axis_Points_2 fed to MeanAxis are: ')
%                    size(Global_Axis_Points_2(:,1,i))
%                     disp('dimensions of Global_Axis_Directions_2 fed to MeanAxis are: ')
%                    size(Global_Axis_Directions_2(:,1,i))
%                end
end

%--------------------------------------------------------------------------

return
end