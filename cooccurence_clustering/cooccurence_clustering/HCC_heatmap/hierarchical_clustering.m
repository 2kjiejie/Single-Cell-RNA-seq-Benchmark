function hierarchical_clustering(data,dim,dist_metric,link_mode,sample_labels)
% hierarchical_clustering(data,dim,dist_metric,link_mode)
%
% parameters:
% 
% dim
%           1 - cluster rows (genes)
%           2 - cluster columns(samples)
% dist_metric
%         'Pearson'     - Pearson correlation
%         'euclidean'   - Euclidean distance (default)
%         'seuclidean'  - Standardized Euclidean distance. Each coordinate
%                         difference between rows in X is scaled by dividing
%                         by the corresponding element of the standard
%                         deviation S=NANSTD(X). To specify another value for
%                         S, use D=PDIST(X,'seuclidean',S).
%         'cityblock'   - City Block distance
%         'minkowski'   - Minkowski distance. The default exponent is 2. To
%                         specify a different exponent, use
%                         D = PDIST(X,'minkowski',P), where the exponent P is
%                         a scalar positive value.
%         'chebychev'   - Chebychev distance (maximum coordinate difference)
%         'mahalanobis' - Mahalanobis distance, using the sample covariance
%                         of X as computed by NANCOV. To compute the distance
%                         with a different covariance, use
%                         D =  PDIST(X,'mahalanobis',C), where the matrix C
%                         is symmetric and positive definite.
%         'cosine'      - One minus the cosine of the included angle
%                         between observations (treated as vectors)
%         'correlation' - One minus the sample linear correlation between
%                         observations (treated as sequences of values).
%         'spearman'    - One minus the sample Spearman's rank correlation
%                         between observations (treated as sequences of values).
%         'hamming'     - Hamming distance, percentage of coordinates
%                         that differ
%         'jaccard'     - One minus the Jaccard coefficient, the
%                         percentage of nonzero coordinates that differ
%         function      - A distance function specified using @, for
%                         example @DISTFUN.
%  
% link_mode:
%        'single'    --- nearest distance (default)
%        'complete'  --- furthest distance
%        'average'   --- unweighted average distance (UPGMA) (also known as
%                        group average)
%        'weighted'  --- weighted average distance (WPGMA)
%        'centroid'  --- unweighted center of mass distance (UPGMC)
%        'median'    --- weighted center of mass distance (WPGMC)
%        'ward'      --- inner squared distance (min variance algorithm)


if ~exist('dim')
    dim=2;
end

if ~exist('dist_metric')
    dist_metric = 'Pearson';
end

if ~exist('link_mode')
    link_mode = 'average';
end



% transpose data according to dim, so that we are always clustering rows
if dim==2
    data = data'; 
end 

if isequal(dist_metric,'Pearson')
    data = data - repmat(mean(data')',1,size(data,2)); % take away mean
    data = data./repmat(std(data')',1,size(data,2));
    dist = data*data'/(size(data,2)-1);
    dist = dist-diag(diag(dist));
    dist = 1-squareform(dist);
else
    dist = pdist(data',dist_metric); 
end

Z = linkage(dist,link_mode); 

if exist('sample_labels')
    dendrogram(Z,size(data,1),'ORIENTATION','left','LABELS',sample_labels)
else
    dendrogram(Z,size(data,1))
end

