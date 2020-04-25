clc
clear all
close all
A = double(imread('x5.jpg'));
A = A / 255;
img_size = size(A);
X = reshape(A, img_size(1) * img_size(2), 3);
fid = fopen('input.txt','w');
for i=1:size(X,1)
    fprintf(fid,'%e\t',X(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
