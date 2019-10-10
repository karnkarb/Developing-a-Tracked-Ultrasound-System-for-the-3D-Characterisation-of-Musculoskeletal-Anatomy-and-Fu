function [slopes, intercepts, fiducial_indices] = get_lines(fiducials)

residual_norm_thresh = 20;

% get all possible combinations of fiducials
combs = combnk(linspace(1, size(fiducials, 1), size(fiducials, 1)), 3);

% test which lines are parallel (((Think it should be collinear here)))
fiducial_indices = zeros(0, 3);
slopes = [];
intercepts = [];

for i = 1:size(combs, 1)
    
    % fit line to the points
    x_coords = [fiducials(combs(i, 1), 1); fiducials(combs(i, 2), 1); fiducials(combs(i, 3), 1)];
    y_coords = [fiducials(combs(i, 1), 2); fiducials(combs(i, 2), 2); fiducials(combs(i, 3), 2)];
    [p, S] = polyfit(x_coords, y_coords, 1); %p is the coefficient
    
    % check if the points are roughly collinear
    if(S.normr < residual_norm_thresh)
        slopes = [slopes, p(1)];
        intercepts = [intercepts, p(2)];
        fiducial_indices = vertcat(fiducial_indices, [combs(i, 1), combs(i, 2), combs(i, 3)]);
    end

end

end

