function [ success_code ] = write_vtk_triamesh_plus_3Ddata(filename, nodes, vector_data, connectivity, id_correction)
%write_vtk_triamesh_plus_3Ddata This writes a vtk triangulated surface plus 3D nodal data
%   v1.0 09.04.2016 (c) Markus Heller m.o.heller@soton.ac.uk
%   
%   This function writes a triangulated surface mesh in vtk format with associated 3D nodal data

success_code = -1; 

%--------------------------------------------------------------------------
% we are given a list of triangles by their 3D coordinates and connectivity
% there is also some data associated to the nodes of the tets
%
% For the connectivity,  the id of the respective nodes might need to 
% be adjusted as vtk indeces start at 0 whilst arrays start with 1 in matlab
% The latter aspect is controlled by the variable id_correction which can
% take values of either 0 or -1 (added to each of the provided connectivity
% data)
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
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

%--------------------------------------------------------------------------
% write the node data
%--------------------------------------------------------------------------

% we need to calculate the number of nodes first
n_nodes = size(nodes,1)

fprintf(fileID, 'POINTS %d float\n', n_nodes);

for i=1:n_nodes;       
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
for i=1:n_trias;       
    fprintf(fileID, '%d ', 3); % a linear tria has 3 nodes
    for j=1:3;
        fprintf(fileID, '%d ', connectivity(i,j)+id_correction);          
    end
    fprintf(fileID, '\n');
end

fprintf(fileID, 'CELL_TYPES %d\n', n_trias); % n_tets
for i=1:n_trias;  
    fprintf(fileID, '5\n'); % a linear tria is VTK cell type no 5
end

%--------------------------------------------------------------------------
% now write the data associated to the nodes
%--------------------------------------------------------------------------
fprintf(fileID, 'POINT_DATA %d\n', n_nodes);     % how many?
fprintf(fileID, 'LOOKUP_TABLE DEFAULT\n');             % how will the data be colored
% fprintf(fileID, 'SCALARS displacements double 3\n');  % name of the stuff
fprintf(fileID, 'VECTORS displacements double\n');  % name of the stuff

for i=1:n_nodes;       
    fprintf(fileID, '%.8f %.8f  %.8f \n', ...
                     vector_data(i,1), ...
                     vector_data(i,2), ...
                     vector_data(i,3)); 
end

%--------------------------------------------------------------------------
% close the file again
%--------------------------------------------------------------------------
fclose(fileID);

success_code = 1;

return;

end