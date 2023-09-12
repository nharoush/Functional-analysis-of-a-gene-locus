function [idHets, id2x]=sortGAP2dosek2(Gap4sort,k)
    %Gap4sort=Kr;
options = statset('UseParallel',1);
[mV,maxID]=max(nansum(Gap4sort,2));
% maxID=find(max(Gap4sort,[],2)==mV,1,'first');

[minV,minID]=min(nansum(Gap4sort,2));

% initialCentroids(1,:)=Gap4sort(minID,:); %check for min and max first, then check other random options
% initialCentroids(2,:)=Gap4sort(maxID,:);
idx=kmeans(Gap4sort,k,'Options',options,'MaxIter',10000,...
    'Display','off','Replicates',30);%, 'Start',initialCentroids);
a=sum(nanmean(Gap4sort(idx==1,:),1));
b=sum(nanmean(Gap4sort(idx==2,:),1));
% c=max(nanmean(Gap4sort(idx==3,:)));
% [m,idNulls]=min([a b c]);
% [m,id2x]=max([a b c]);
[m,idHets]=min([a b]);
[m,id2x]=max([a b]);

idHets=find(idx==idHets );
% idNulls=(find(idx==idNulls));
id2x=(find(idx==id2x));

end
