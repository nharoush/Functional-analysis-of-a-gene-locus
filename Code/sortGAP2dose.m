function [idHets, idNulls,id2x]=sortGAP2dose(Gap4sort)
%This is specifically when k=3; instead of default centroid initialization, use our knowledge
%about the 2x-1x clusters and initialize to the largest and smallest max
%for each  time sample

k=3;
options = statset('UseParallel',1);
% mV=max(Gap4sort, [],'all');
% maxID=find(max(Gap4sort,[],2)==mV,1,'first');
% 
% minV=min(min(Gap4sort));
% minID=find(min(Gap4sort,[],2)==minV,1,'first');
% initialCentroids(1,:)=Gap4sort(minID,:);
% initialCentroids(3,:)=Gap4sort(maxID,:);
% initialCentroids(2,:)=mean([Gap4sort(maxID,:),Gap4sort(minID,:)]);
% 
% initialCentroids=repmat(initialCentroids, 1,1,30);
% noise=0.1*randn(size(initialCentroids));
% initialCentroids=initialCentroids+noise;
    idx=kmeans(Gap4sort,k,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',30);%, 'Start',initialCentroids);%,'Distance', 'Cosine');

a=max(nanmean(Gap4sort(idx==1,:)));
b=max(nanmean(Gap4sort(idx==2,:)));
c=max(nanmean(Gap4sort(idx==3,:)));
[m,idNulls]=min([a b c]);
[m,id2x]=max([a b c]);

idHets=find(idx~=idNulls & idx~=id2x);
idNulls=(find(idx==idNulls));
id2x=(find(idx==id2x));
end
