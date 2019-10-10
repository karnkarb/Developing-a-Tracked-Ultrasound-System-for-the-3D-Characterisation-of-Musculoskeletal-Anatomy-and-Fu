function [ M ] = create_4by4_tmatrix ( R, ori )

% the function establishes a 4x4 homogeneous transformation matrix M.
% Inputs a a 3x3 rotation matrix R and the location of the origin of the 
% local coordinate system, ori.
% M is supposed to be applied to COLUMN vectors, and the rotation matrix R 
% is supposed to be arranged accordingly (i.e. direction vectors are to be 
% provided in ROWs.
% The origin, ori, can be provided as either row or column vector and will 
% be transposed if necessary 

% CHECKS
[ n_rows, n_cols ] = size(R);
if ( n_rows ~= 3 || n_cols ~= 3 )
    msg_str = sprintf('ERROR in create_4by4_tmatrix: R is not 3x3 but of dims: %d x %d\n', ...
        n_rows, n_cols);
    R
    disp(msg_str);
    M=[];
    return;
end    

% initialise M as 4x4 I matrix
M = eye(4);

% position the rotation sub-matrix R
M(1:3,1:3) = R;

% do we need to transpose ori?
if ( size(ori,1) == 3 && size(ori,2) == 1 )
    ori = ori';
end

% prepare the translation to the origin of the local CS
t = -ori';
t = R*t;

% compose the full 4x4 matrix from R and t
M(1:3,4) = t;

end