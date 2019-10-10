function [img_output] = noise_cleaning_linear(img_input, k, k2)

se = strel('rectangle', k);
se2 = strel('disk', k2);
se3 = strel ('line', 12,2);

img_output = imopen(img_input, se);
img_output = imclose(img_output, se);

img_output = imopen(img_output, se2);
img_output = imclose(img_output, se2);

img_output = imdilate(img_output,se3);

end

 