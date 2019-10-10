function [ success_code ] = write_vtk_points(filename, points)
%write_vtk_points This writes a set of vtk points
%   v1.0 29.03.2015 (c) Markus Heller m.o.heller@soton.ac.uk
%   
%   This function writes vtk point data

success_code = -1; 

% DEBUG: show what we have got
% size(points)
% size(connectivity)

%--------------------------------------------------------------------------
% open the file for writing
%--------------------------------------------------------------------------
fileID = fopen(filename,'w');

%--------------------------------------------------------------------------
% write the header info
%--------------------------------------------------------------------------
fprintf(fileID, '# vtk DataFile Version 3.0\n');
fprintf(fileID, '%s\n', filename);
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET POLYDATA\n');

%--------------------------------------------------------------------------
% write the point data
%--------------------------------------------------------------------------
fprintf(fileID, 'POINTS %d float\n', size(points, 1));

for i = 1:size(points, 1);    
    fprintf(fileID, '%.8f %.8f %.8f\n', ...
                     points(i,1),  points(i,2),  points(i,3));
end

%--------------------------------------------------------------------------
% close the file again
%--------------------------------------------------------------------------
fclose(fileID);

success_code = 1;

return;

end