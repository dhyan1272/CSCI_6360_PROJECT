function idx = findClosestCentroids(X, centroids)

idx = zeros(size(X,1), 1);
ind = zeros(1,1);
dist=zeros(size(centroids,1),1);
m=size(X,1);
n=size(centroids,1);
for i=1:m
    for j=1:n    
        dist(j)=sum((X(i,:)-centroids(j,:)).^2);
        ind=find(dist==min(dist));
    end
    idx(i,1)=ind(1,1);
    
end   
end

