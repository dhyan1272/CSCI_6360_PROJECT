clc
clear all
close all
A = double(imread('x2.jpg'));
A = A / 255; % Divide by 255 so that all values are in the range 0 - 1
img_size = size(A);
X = reshape(A, img_size(1) * img_size(2), 3);
K = 3; 
max_iters = 1;
initial_centroids = kMeansInitCentroids(X, K);
centroids=initial_centroids;
for i=1:max_iters
   
    idx = findClosestCentroids(X, centroids);
    centroids = computeCentroids(X, idx, K);
end

fprintf('\nApplying K-Means to compress an image.\n\n')
idx = findClosestCentroids(X, centroids);
X_recovered = centroids(idx,:);

% Reshape the recovered image into proper dimensions
X_recovered = reshape(X_recovered, img_size(1), img_size(2), 3);

% Display the original image 
subplot(1, 2, 1);
imagesc(A); 
title('Original');

% Display compressed image side by side
subplot(1, 2, 2);
imagesc(X_recovered)
title(sprintf('Compressed, with %d colors.', K));

function centroids = kMeansInitCentroids(X, K)
centroids = zeros(K, size(X, 2));
randidx = randperm(size(X, 1));
centroids = X(randidx(1:K), :);
end
