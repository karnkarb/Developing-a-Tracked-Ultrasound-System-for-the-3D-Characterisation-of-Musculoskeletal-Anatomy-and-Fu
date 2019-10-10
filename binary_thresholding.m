function output_img = binary_thresholding(input_img, t)

%find optimal thresholding
thre = (multithresh(input_img) - t)/255.0;

%normalise the image
I = input_img/255.0;

%binarize the image
output_img = im2bw(I, thre);


end

