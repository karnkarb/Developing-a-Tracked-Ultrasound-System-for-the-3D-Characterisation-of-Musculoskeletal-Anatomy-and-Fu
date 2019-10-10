function slope_indices = get_parallel_three_lines(slopes, fiducial_indices)

min_slope_difference = 0.1;
slope_indices = zeros(1, 3);

for i = 1:size(slopes, 2) - 2
    for j = (i+1):size(slopes, 2) -1
        for k = (i+2):size(slopes, 2)
        % check if fiducials are all different
        if isempty ([intersect(fiducial_indices(i, :), fiducial_indices(j, :)) ...
            intersect(fiducial_indices(i, :), fiducial_indices(k, :))...
            intersect(fiducial_indices(j, :), fiducial_indices(k, :))])
        diff_slope1 = abs(slopes(i)-slopes(k));
        diff_slope2 = abs(slopes(i)-slopes(j));
        diff_slope3 = abs(slopes(j)-slopes(k));
     % check the difference between slopes
        if (diff_slope1<min_slope_difference)&&(diff_slope2<min_slope_difference)&&(diff_slope3<min_slope_difference)
            slope_indices = [i, j, k];
            return
        else 
            slope_indices = [ 0 0 0];
        end
        else
            continue
        end
        
        end
                
    end
end
end