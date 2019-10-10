function [img_output] = noise_cleaning(img_input, k, k2)

se = strel('rectangle', k);
se2 = strel('disk', k2);
se3 = strel('line',5,0);

img_output = imopen(img_input, se);
img_output = imclose(img_output, se);
img_output = imdilate(img_output,se3);

img_output = imopen(img_output, se2);
img_output = imclose(img_output, se2);

end

