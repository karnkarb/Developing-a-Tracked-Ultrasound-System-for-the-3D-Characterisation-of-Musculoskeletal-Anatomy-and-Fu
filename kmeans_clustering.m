function [img_output, img_output2, colours_to_remove] = kmeans_clustering(img_input, k)

%create feature matrix
[r, c] = find(img_input==1);
X = horzcat(r, c);

%apply k means
idx = kmeans(X, k);
img_output = idx;

%generate k random colours
colours = rand(k, 3);

%create coloured image
img_output = zeros(size(img_input, 1), size(img_input, 2), 3);
for i = 1:size(X, 1)
   img_output(X(i, 1), X(i, 2), :) = colours(idx(i), :);
end

%remove connected components
colours_to_remove = zeros(0, 3);
R_margin = size(img_input, 1);
C_margin = size(img_input, 2);

for i = 1:size(X, 1)
    pixel_colour = img_output(X(i, 1), X(i, 2), :);
    pixel_colour = reshape(pixel_colour, 1,3);
    % check if already removed
    if all(ismember(pixel_colour, colours_to_remove))
        continue
    end
    
    %check if 2 different components are connected
    coordinates = X(i, :);
    x = coordinates(1);
    y = coordinates(2);
    
    for r = max(1, x-1):min(R_margin, x+1)
        for c = max(1, y-1):min(C_margin, y+1)
            neighbour_colour = img_output(r, c, :);
            neighbour_colour = reshape(neighbour_colour, 1, 3);
            if (all(neighbour_colour == pixel_colour) || sum(neighbour_colour) == 0)
                continue
            end
            colours_to_remove = vertcat(colours_to_remove, neighbour_colour);
            colours_to_remove = vertcat(colours_to_remove, pixel_colour);
        end
    end
        
end

img_output2 = img_output;

for r = 1:size(img_output2, 1)
    for c = 1:size(img_output2, 2)
        if all(ismember(img_output2(r, c, :), colours_to_remove))
            img_output2(r, c, :) = [0, 0, 0];
        end
    end
end



end

