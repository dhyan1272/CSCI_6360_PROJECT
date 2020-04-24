function centroids = computeCentroids(X, idx, K)
[m n] = size(X);
centroids = zeros(K, n);
for i=1:K
    ind=find(idx==i);
    n_cen(i,:)=(sum(X(ind,:)))/size(ind,1);
end
centroids=n_cen;
end

