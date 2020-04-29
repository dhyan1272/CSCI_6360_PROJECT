clc
clear all
close all

%Original Image
A = double(imread('x3.jpg'));
A = A / 255; % Divide by 255 so that all values are in the range 0 - 1
img_size = size(A);

%Image from C program
fileID=fopen('output.bin');
X_bin=fread(fileID,'double');
fclose(fileID);
for i=1:img_size(1) * img_size(2) 
    X1(i,:)=X_bin((i-1)*3+1:3*i);
end
A2 = reshape(X1, img_size(1), img_size(2), 3);

%Display original image
subplot(1, 2, 1);
imagesc(A); 
title('Original');

% Display compressed image
subplot(1, 2, 2);
imagesc(A2)
title(sprintf('Compressed image'));
