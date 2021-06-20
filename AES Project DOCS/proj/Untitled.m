rgbImage=imread('ex.jpeg');
red = rgbImage(:,:,1);
allBlack = zeros(size(rgbImage, 1), size(rgbImage, 2), 'uint8');
just_red = cat(3, red, allBlack, allBlack);
imshow(just_red)
just_red
%imshow(i);