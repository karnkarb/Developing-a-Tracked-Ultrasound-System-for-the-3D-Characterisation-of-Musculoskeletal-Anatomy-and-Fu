function [ transformed_data ] = apply_4by4_tmatrix ( M, data, varargin )

% This function applies a 4x4 homogenous transformation matrix M 
% to an mx3 or 3xm data matrix, data.
% M is supposed to be applied to data presented as COLUMN vectors,
% but the data matrix can be arranged in either column or row format.
%
% NOTE:
% =====
% We could not possibly figure out what to do with a 3x3 matrix of points
%  
% --> in that case (or otherwise) we can provide a 3rd input argument:
%
%     varargin == 1: COLUMN VECTOR MODE
%              == 2: ROW VECTOR MODE
%--------------------------------------------------------------------------
% v1.1 21.10.2018 (c) Markus Heller m.o.heller@soton.ac.uk
%--------------------------------------------------------------------------

[ n_rows, n_columns ] = size(data);

% first: check which of the rows/columns is of size 3; 
%        if both rows and columns are of size 3, assume row data, unless 
%        instructed otherwise
if ( n_rows ~= 3 && n_columns ~= 3 )
    msg_str = sprintf('ERROR: expect either n_rows or n_columns to be 3 - but we havd %d rows and %d columns!!!\nWill copy the data, not transform it!!!', ...
                       n_rows, n_columns );
    disp(msg_str);
    
    transformed_data = data;
    return;
end

if ( n_rows == 3 && n_columns == 3 )
    msg_str = sprintf('WARNING: 3 by 3 matrix!\nCannot uniquely identify vector mode\nWill assume ROW vectors - check if that is indeed correct!!!');
    disp(msg_str);
    
    if ( nargin == 3 )
        if ( varargin{1} == 1 ) % COLUMN VECTOR mode
            row_input_data = 0;
            msg_str = sprintf('pooh: 3rd argument used to enter COLUMN VECTOR mode !!!');
            disp(msg_str);
        else
            row_input_data = 1;
            data = data'; 
            msg_str = sprintf('pooh: 3rd argument used to enter ROW VECTOR mode !!!');
            disp(msg_str);
        end
    else 
        msg_str = sprintf('WARNING: 3 by 3 matrix!\nCannot uniquely identify vector mode\nWill assume ROW vectors - check if that is indeed correct!!!');
        disp(msg_str);

        row_input_data = 1;  
        data = data'; 
    end
else
    if ( n_columns == 3 )
        row_input_data = 1;   
        data = data'; 
    else
        row_input_data = 0;
    end
end
    
n_pts = size(data,2);
% we are now sure that we have a data matrix 3xn_pts ----------------------

% create a homogenous version of the data (i.e. 4xm)
hom_data = [ data; ones(1, n_pts) ];
    
% now transform the data
transformed_data = M * hom_data;

% get rid of the 4th row
transformed_data(4,:)=[];

% if we had row vectors as input, we need to transform the results to 
% match that data representation
if ( row_input_data )
    transformed_data = transformed_data';
end

end