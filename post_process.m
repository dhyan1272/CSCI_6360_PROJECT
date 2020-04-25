clc
clear all
close all

A = double(imread('x5.jpg'));
A = A / 255; % Divide by 255 so that all values are in the range 0 - 1
img_size = size(A);
X = reshape(A, img_size(1) * img_size(2), 3);


A2= importdata('output.txt');
A2 = reshape(A2, img_size(1), img_size(2), 3);
imagesc(A); 

subplot(1, 2, 1);
imagesc(A); 
title('Original');

% Display compressed image side by side
subplot(1, 2, 2);
imagesc(A2)
title(sprintf('Compressed image'));
