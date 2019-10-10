function [ success_code ] = write_vtk_multipolyline( filename, polylines )
%write_vtk_multipolyline This writes a collection of vtk polylines
%   v1.0 10.04.2015 (c) Markus Heller m.o.heller@soton.ac.uk
%   
%   This function writes data for a collection of vtk polylines
%   Polylines can either be open, or closed.
%
% inputs:
% =======
% filename      string, filename
% polylines     array of structures with fields:
%               polylines.vertices  -- n_vertices x 3
%               polylines.is_closed -- scalar value: 0 if open, 1 if closed          
%--------------------------------------------------------------------------

DEBUG = 1;

success_code = -1; 

n_polylines = size(polylines,2);
if DEBUG
   info_str = sprintf('Will write %d vtk polylines',  n_polylines);
   disp(info_str);
end

%--------------------------------------------------------------------------
% we are given a list of polylines by their 3D coordinates
% connectivity is simple, just the id of the respective point minus 1(-1)
% (vtk indices start at 0 whilst arrays start with 1 in matlab)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% open the file for writing
%--------------------------------------------------------------------------
fileID = fopen(filename,'w');

%--------------------------------------------------------------------------
% write the header info
%--------------------------------------------------------------------------
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, '%s\n', filename);
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET POLYDATA\n');

%--------------------------------------------------------------------------
% write the point data
%--------------------------------------------------------------------------

% we need to calculate the number of points first
n_points = 0;
for i=1:n_polylines
    if ( polylines(i).is_closed )
        n_points = n_points + size(polylines(i).vertices,1) - 1;
    else
        n_points = n_points + size(polylines(i).vertices,1);
    end
end
    
fprintf(fileID, 'POINTS %d float\n', n_points);

n_vertices = 0;
for i=1:n_polylines
    n_vertices = n_vertices + size(polylines(i).vertices,1);
    if ( polylines(i).is_closed )
        for j=1:size(polylines(i).vertices,1)-1
            fprintf(fileID, '%.8f %.8f %.8f\n', ...
                        polylines(i).vertices(j,1), ... 
                        polylines(i).vertices(j,2), ... 
                        polylines(i).vertices(j,3) );
        end
    else
        for j=1:size(polylines(i).vertices,1)
            fprintf(fileID, '%.8f %.8f %.8f\n', ...
                        polylines(i).vertices(j,1), ... 
                        polylines(i).vertices(j,2), ... 
                        polylines(i).vertices(j,3) );
        end        
    end
end

%--------------------------------------------------------------------------
% write the connectivity info for the line - how to connect the points ...
%--------------------------------------------------------------------------
vtksize = n_polylines + n_vertices;
fprintf(fileID, 'LINES %d %d\n', n_polylines, vtksize);

p_index = 0;
% loop over all lines ...
for i=1:n_polylines  
    if ( polylines(i).is_closed )
        fprintf(fileID, '%d ', size(polylines(i).vertices,1));
        for j=1:size(polylines(i).vertices,1)-1
            fprintf(fileID, '%d ', p_index );   
            p_index = p_index + 1;
        end
        fprintf(fileID, '%d ', p_index-(size(polylines(i).vertices,1)-1) ); 
        fprintf(fileID, '\n');
    else
        fprintf(fileID, '%d ', size(polylines(i).vertices,1));
        for j=1:size(polylines(i).vertices,1)
            fprintf(fileID, '%d ', p_index );   
            p_index = p_index + 1;
        end
        fprintf(fileID, '\n');
    end
end

%--------------------------------------------------------------------------
% close the file again
%--------------------------------------------------------------------------
fclose(fileID);

success_code = 1; % all good - set success code to 1

return;

end