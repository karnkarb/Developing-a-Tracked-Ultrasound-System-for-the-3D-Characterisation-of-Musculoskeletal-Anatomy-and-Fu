function [coordinates] = get_coordinates(best_lines, fiducial_indices, centroids)

%create empty matric for coordinates to go in
coordinates = [];

%find the two best lines output
line1 = best_lines(1);
line2 = best_lines(2);
line3 = best_lines(3);

%find the three fiducial indeces each line is associated with
line1 = fiducial_indices(line1,:);
line2 = fiducial_indices(line2,:);
line3 = fiducial_indices(line3,:);

%the number of each fiducial
f1 = line1(1);
f2 = line1(2);
f3 = line1(3);
f4 = line2(1);
f5 = line2(2);
f6 = line2(3);
f7 = line3(1);
f8 = line3(2);
f9 = line3(3);


coordinate_1 = centroids(f1,:);
coordinate_2 = centroids(f2,:);
coordinate_3 = centroids(f3,:);
coordinate_4 = centroids(f4,:);
coordinate_5 = centroids(f5,:);
coordinate_6 = centroids(f6,:);
coordinate_7 = centroids(f7,:);
coordinate_8 = centroids(f8,:);
coordinate_9 = centroids(f9,:);

coordinates(1,:) = coordinate_1;
coordinates(2,:) = coordinate_2;
coordinates(3,:) = coordinate_3;
coordinates(4,:) = coordinate_4;
coordinates(5,:) = coordinate_5;
coordinates(6,:) = coordinate_6;
coordinates(7,:) = coordinate_7;
coordinates(8,:) = coordinate_8;
coordinates(9,:) = coordinate_9;



end
















