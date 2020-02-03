function [perm_r,perm_c] = HCC_heatmap(data,mode_of_sim)


% This part is for visualization the training data using hierarchical clustering provided by the matlab statistical toolbox
if exist('mode_of_sim')
    if mode_of_sim=='r'
        eucD = pdist(data,'euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
        perm_c = 1:size(data,2);
    elseif mode_of_sim=='c'
        eucD = pdist(data','euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
        perm_r = 1:size(data,1);
    else
        eucD = pdist(data,'euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
        eucD = pdist(data','euclidean');
        clustTreeEuc = linkage(eucD,'average');
        [H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
    end
    imagesc(data(perm_r,perm_c));
    return
end

eucD = pdist(data,'euclidean');
clustTreeEuc = linkage(eucD,'average');
[H,T,perm_r] = dendrogram(clustTreeEuc,size(data,1));
eucD = pdist(data','euclidean');
clustTreeEuc = linkage(eucD,'average');
[H,T,perm_c] = dendrogram(clustTreeEuc,size(data,2));
imagesc(data(perm_r,perm_c));
