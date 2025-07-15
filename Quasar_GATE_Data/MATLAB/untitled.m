id = fopen('Phase70_cropped.raw');
img = fread(id,'short');
imgsize = size(img);
rows = 237;
clos = 158;
nums = imgsize(1)/rows/clos; 
img = reshape(img,[rows,clos,nums]);
single_image = reshape(img(:,:,3),[rows,clos]);
single_show = uint8(single_image);
imshow(single_show)