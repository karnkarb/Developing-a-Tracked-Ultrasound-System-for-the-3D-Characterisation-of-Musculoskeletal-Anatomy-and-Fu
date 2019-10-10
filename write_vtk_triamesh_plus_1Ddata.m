function [ success_code ] = write_vtk_triamesh_plus_1Ddata(filename, nodes, scalar_data, scalar_data_name, connectivity, id_correction)
%write_vtk_triamesh_plus_1Ddata This writes a vtk tria surface plus data
%   v1.0 26.03.2016 (c) Markus Heller m.o.heller@soton.ac.uk
%   v1.1 16.02.2019 (c) Markus Heller m.o.heller@soton.ac.uk
%                       ENH more robust handling of the content of scalar_data_name
%   This function writes a tria surface mesh with associated scalar data in vtk format
%--------------------------------------------------------------------------
% we are given a list of triangles by their 3D coordinates and connectivity
% there is also some data associated to the nodes of the tets
%
% For the connectivity,  the id of the respective nodes might need to 
% be adjusted as vtk indeces start at 0 whilst arrays start with 1 in matlab
% The latter aspect is controlled by the variable id_correction which can
% take values of either 0 or -1 (added to each of the provided connectivity
% data)
%
% inputs:
% =======
% filename          string providing the filename (incl. full path)
% nodes             n_nodes x 3 array with 3D node coordinates
% scalar_data       vector with scalar data associated to the nodes
% scalar_data_name  brief string providing name for scalar data;
%                   should not containy any white space characters, but 
%                   if there are any, they will be removed
% connectivity      n_trias x 3 array with triangle connectivity info
% id_correction     0 if the connectivity is specified with 0 indexed nodes
%                   (i.e. node indices start from 0),
%                   -1 if the connectivity is specified with 1 indexed nodes
%                   (i.e. node indices start from 1)
%--------------------------------------------------------------------------

success_code = -1; 

%--------------------------------------------------------------------------
% open the file for writing
%--------------------------------------------------------------------------
fileID = fopen(filename,'w');

if ( fileID == -1 )
	error('could not open file: %s for writing!', filename);
end

%--------------------------------------------------------------------------
% write the header info
%--------------------------------------------------------------------------
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, '%s\n', filename);
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

%--------------------------------------------------------------------------
% write the node data
%--------------------------------------------------------------------------

% we need to calculate the number of nodes first
n_nodes = size(nodes,1)

fprintf(fileID, 'POINTS %d float\n', n_nodes);

for i=1:n_nodes      
    fprintf(fileID, '%.8f %.8f %.8f\n', ...
                     nodes(i,1), ... 
                     nodes(i,2), ... 
                     nodes(i,3) ); 
end

%--------------------------------------------------------------------------
% write the connectivity info for the tets - how to connect the nodes ...
%--------------------------------------------------------------------------
n_trias = size(connectivity,1)  
vtksize = 4 * n_trias;
fprintf(fileID, 'CELLS %d %d\n', n_trias, vtksize);

% loop over all tets ...
for i=1:n_trias     
    fprintf(fileID, '%d ', 3); % a linear tria has 3 nodes
    for j=1:3
        fprintf(fileID, '%d ', connectivity(i,j)+id_correction);          
    end
    fprintf(fileID, '\n');
end

fprintf(fileID, 'CELL_TYPES %d\n', n_trias); % n_tets
for i=1:n_trias  
    fprintf(fileID, '5\n'); % a linear tria is VTK cell type no 5
end

%--------------------------------------------------------------------------
% now write the data associated to the nodes
%--------------------------------------------------------------------------
fprintf(fileID, 'POINT_DATA %d\n', n_nodes);     % how many?

% make sure the name of the scalars does not cause trouble because of " "
scalar_data_name = regexprep(scalar_data_name,'\s*','');

fprintf(fileID, 'SCALARS %s double 1\n', scalar_data_name);  % name of the scalar data
fprintf(fileID, 'LOOKUP_TABLE default\n');       % how will the data be colored

for i=1:n_nodes       
    fprintf(fileID, '%.8f\n', scalar_data(i,1)); 
end

%--------------------------------------------------------------------------
% close the file again
%--------------------------------------------------------------------------
fclose(fileID);

success_code = 1;

return

end