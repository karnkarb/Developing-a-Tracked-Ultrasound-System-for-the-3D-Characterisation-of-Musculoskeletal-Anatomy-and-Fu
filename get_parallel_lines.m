function slope_indices = get_parallel_lines(slopes, fiducial_indices)

min_slope_difference = 0.1;
slope_indices = zeros(1, 2);

for i = 1:size(slopes, 2) - 1
    for j = (i+1):size(slopes, 2)
        
        % check if fiducials are all different
        if(intersect(fiducial_indices(i, :), fiducial_indices(j, :)))
            continue
        end
        
        % check the difference between slopes
        if(abs(slopes(i) - slopes(j)) < min_slope_difference)
            min_slope_difference = abs(slopes(i) - slopes(j));
            
            slope_indices = [i, j];
        end
        
    end
end

end