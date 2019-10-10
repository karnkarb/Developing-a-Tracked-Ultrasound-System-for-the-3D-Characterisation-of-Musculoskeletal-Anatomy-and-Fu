function fiducial_centroids = get_fiducial_centroids(input_img)

% binarize image
epsilon = 1e-3;
binary_img = im2bw(input_img, epsilon);

% find connected components
connected_comps = bwconncomp(binary_img);
comp_centroids = regionprops(connected_comps, 'Centroid');

% get centroids as coordinates matrix
centroids_cell = struct2cell(comp_centroids);
fiducial_centroids = reshape(cell2mat(centroids_cell), 2, size(centroids_cell, 2))';

end

