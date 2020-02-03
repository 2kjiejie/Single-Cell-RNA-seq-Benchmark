function [cell_labels, gene_labels, cell_cluster_best_snr, signature_for_which_celltype] = cooccurance_clustering(data, min_expressed_cells, min_pathway_size, min_population_size, snr_merge_threshold)

% data = double(counts_data~=0); % data here is the binary version
fprintf(1, 'Processing data subset with %d genes and %d cells: \n', size(data,1), size(data,2));
gene_labels = zeros(size(data,1),1);
cell_labels = zeros(size(data,2),1);
cell_cluster_best_snr = [];
signature_for_which_celltype = [];

% filter out genes detected in too few cells
tic
keep_ind = find(sum(data~=0,2)>=min_expressed_cells & sum(data==0,2)>=min_expressed_cells);  % only keep genes that express in >=10 cells
data = data(keep_ind,:);
fprintf(1, 'Remove genes detected in <%d cells. Remaining %d genes.  ', min_expressed_cells, size(data,1));
toc
if length(keep_ind)==0
    return
end


% Finding a threshold for the gene-gene chi2 values using random permute
tic
fprintf(1, 'Iterate 10 random permutations for gene-gene similarity threshold ... %4d ', 0)
rand_iter = 10;
tmp = [];
for i=1:rand_iter
    fprintf(1, '\b\b\b\b\b%4d ', i)
    perm_data = zeros(size(data));
    [~,I] = sort(rand(size(data)),2);
    perm_data = data((I-1)*size(data,1)+repmat((1:size(data,1))',1,size(data,2)));
    clear I
    A_plus_B_vector = sum(perm_data~=0,2);  
    C_plus_D_vector = sum(perm_data==0,2); 
    A_plus_C_vector = sum(perm_data~=0,2)';
    B_plus_D_vector = sum(perm_data==0,2)';
    chi2_denominator = ((A_plus_B_vector*A_plus_C_vector).*(C_plus_D_vector*B_plus_D_vector));
    A = perm_data*perm_data';
    B = A_plus_B_vector - A; % B = A_plus_B - A;  % B2 = double(data~=0)*double(data==0)';
    C = A_plus_C_vector - A; % C = A_plus_C - A;  % C2 = double(data==0)*double(data~=0)';
    D = C_plus_D_vector - C; % D = C_plus_D - C;  % D2 = double(data==0)*double(data==0)';
    clear A_plus_B_vector C_plus_D_vector A_plus_C_vector B_plus_D_vector
    AD = A.*D; 
    clear A D
    BC = B.*C;
    clear B C
    AD_BC = AD-BC;
    clear AD BC
    chi2_numerator = sign(AD_BC).*(AD_BC.^2);
    clear AD_BC
    rand_chi2_values = chi2_numerator./chi2_denominator; 
    clear chi2_numerator chi2_denominator
    rand_chi2_values = rand_chi2_values - diag(diag(rand_chi2_values));                                         % get rid of diagnal which should be all 1s
    rand_chi2_values(isnan(rand_chi2_values))=0;                                                                % get rid of NaN
    rand_chi2_values(rand_chi2_values<0)=0;                                                                     % zero out the negative value
    rand_chi2_values(rand_chi2_values==1) = 0;                                                                  % reduce 1s to the next largest value
    s = squareform(rand_chi2_values);
    clear rand_chi2_values
    Y = sort(s,'descend');
    clear s
    if length(Y)==1
        tmp = [tmp,Y(1)];  % this will only happen if there are just two genes
    else
        tmp = [tmp,Y(2)];  % pick the second largest rand chi2 value
    end
    clear Y
end
best_random_chi2_value = median(tmp);
toc
clear perm_data





% Gene-gene similarity using all cells
tic
fprintf(1, 'Compute gene-gene similarity ... ')
num_cells = size(data,2); 
num_genes = size(data,1);
A_plus_B_vector = sum(data~=0,2);  % A_plus_B = repmat(sum(data~=0,2),1,size(data,1));  
C_plus_D_vector = sum(data==0,2);  % C_plus_D = repmat(sum(data==0,2),1,size(data,1));
A_plus_C_vector = sum(data~=0,2)'; % A_plus_C = repmat(sum(data~=0,2)',size(data,1),1);
B_plus_D_vector = sum(data==0,2)'; % B_plus_D = repmat(sum(data==0,2)',size(data,1),1);
chi2_denominator = ((A_plus_B_vector*A_plus_C_vector).*(C_plus_D_vector*B_plus_D_vector));
A = data*data';
B = A_plus_B_vector - A; % B = A_plus_B - A;  % B2 = double(data~=0)*double(data==0)';
C = A_plus_C_vector - A; % C = A_plus_C - A;  % C2 = double(data==0)*double(data~=0)';
D = C_plus_D_vector - C; % D = C_plus_D - C;  % D2 = double(data==0)*double(data==0)';
clear A_plus_B_vector C_plus_D_vector A_plus_C_vector B_plus_D_vector
AD = A.*D; 
clear A D
BC = B.*C;
clear B C
AD_BC = AD-BC;
clear AD BC
chi2_numerator = sign(AD_BC).*(AD_BC.^2);
clear AD_BC
all_chi2_values = chi2_numerator./chi2_denominator; 
clear chi2_numerator chi2_denominator
all_chi2_values = all_chi2_values - diag(diag(all_chi2_values));                                        % get rid of diagnal which should be all 1s
all_chi2_values(isnan(all_chi2_values))=0;                                                              % get rid of NaN
all_chi2_values(all_chi2_values<0)=0;                                                                   % zero out the negative value
all_chi2_values(all_chi2_values==1) = max(all_chi2_values(all_chi2_values~=1));                         % reduce 1s to the next largest value
toc


% defind the threshold
[Y,I] = sort(all_chi2_values,2,'descend');
mu = mean(Y(Y(:,1)<best_random_chi2_value,1)); 
sigma = std(Y(Y(:,1)<best_random_chi2_value,1));
threshold = mu + sigma;
if isnan(threshold) || threshold<=0 || isinf(threshold)
    threshold = min(Y(:,1));
end
clear Y I mu sigma best_random_chi2_value

% Use the threshold to select genes share co-occurrence patterns/clusters
tic
fprintf(1, 'Create gene-gene graph for clustering genes ... \n');

    % lower_triangle_adj = tril(all_chi2_values,-1);
    % lower_triangle_adj(lower_triangle_adj<threshold)=0;
    % Graph2Binary(lower_triangle_adj,'G')
    % niter = 20;
    % [c,Q,labels,communities] = LouvainfromBin('G.bin',niter);
    % clear lower_triangle_adj all_chi2_values
    
    A = all_chi2_values;
    A(A<threshold) = 0;
    if sum(sum(A~=0))==0
        labels = zeros(1,size(A,1));
    else
        graph_file = write_symmetric_adj_into_graph_file(A, 'G.txt');
        [labels] = clustering_graph_by_modularity(graph_file);
    end
    clear all_chi2_values A
    
[N] = hist(labels, 1:max(labels)); 
gene_cluster_ID = find(N>=min_pathway_size); 
selected_genes = find(ismember(labels,gene_cluster_ID));
fprintf(1, 'Gene-gene graph contains %d pathways, %d genes in total\n', length(unique(labels(selected_genes))), length(selected_genes));

if isempty(selected_genes)
    toc
    return
end
selected_gene_labels = labels(selected_genes);
gene_labels(keep_ind(selected_genes)) = selected_gene_labels;
gene_labels = standardize_idx(gene_labels);
toc


% Use the selected genes to construct cell-cell graph and clusters
fprintf(1, 'Create cell-cell graph for clustering cells ... \n');
selected_data = data(selected_genes, :);
selected_data_gene_cluster_mean = grpstats(selected_data, selected_gene_labels);
if size(data,2)>=60
    k=30;
elseif size(data,2)>=10
    k = round(size(data,2)/2);
elseif k>=6
    k=5;
else
    k = size(data,2)-1;
end

distance = 'euclidean';
if size(selected_data_gene_cluster_mean,1)>=3
    [IDX,D] = knnsearch(selected_data_gene_cluster_mean',selected_data_gene_cluster_mean','k',k+1,'Distance',distance);
    IDX(:,1) = []; D(:,1) = [];
else
    tmp = [];
    for ss=1:max(selected_gene_labels)
        ind_tmp = find(selected_gene_labels==ss);
        ind_tmp1 = ind_tmp(1:round(length(ind_tmp)/2));
        ind_tmp2 = ind_tmp(round(length(ind_tmp)/2)+1:end);
        tmp = [tmp; mean(selected_data(ind_tmp1,:),1);mean(selected_data(ind_tmp2,:),1)];
    end
    [IDX,D] = knnsearch(tmp',tmp','k',k+1,'Distance',distance);
    IDX(:,1) = []; D(:,1) = [];
end
G = knn2jaccard(IDX);

    % Graph2Binary(G,'G')
    % niter = 20;
    % [c,Q,labels,communities] = LouvainfromBin('G.bin',niter);
    
    A = G+G';
    graph_file = write_symmetric_adj_into_graph_file(A, 'G.txt');
    [labels] = clustering_graph_by_modularity(graph_file);    

cell_labels = labels(:);
fprintf(1, 'Cell-cell graph contains %d cell types by community detection\n', length(unique(cell_labels)));
toc


% decide whether these clusters are good, whether we should merge
if max(cell_labels)==1
    return
end

% merge tiny cell clusters with the non-tiny cluster they are most connected to
cluster_size = hist(cell_labels,0:max(cell_labels)); 
cluster_size(1)=[];
tiny_clusters = find(cluster_size<min_population_size);
tiny_clusters_merge_to = zeros(size(tiny_clusters));
for i=1:length(tiny_clusters)
    tmp = sum(grpstats(full(A(:,find(cell_labels==tiny_clusters(i)))),cell_labels,'max')',1);
    tmp(tiny_clusters) = 0;
    if max(tmp)==0
        pairwise_cluster_center_dist = squareform(pdist(grpstats(selected_data_gene_cluster_mean',cell_labels)));
        tmp = pairwise_cluster_center_dist(tiny_clusters(i),:);
        tmp(tiny_clusters) = Inf;
        [~,tiny_clusters_merge_to(i)] = min(tmp);
    else
        [~,tiny_clusters_merge_to(i)] = max(tmp);
    end
end
for i=1:length(tiny_clusters)
    cell_labels(cell_labels==tiny_clusters(i)) = tiny_clusters_merge_to(i);
end
fprintf(1, 'Cell-cell graph contains %d cell types after merging tiny cell clusters\n', length(unique(cell_labels)));

while 1
    tmp = [];
    for i=1:max(cell_labels)
        for j=i+1:max(cell_labels)
            best_snr  = max(abs(get_correlations(selected_data_gene_cluster_mean(:,ismember(cell_labels,[i,j])), cell_labels(ismember(cell_labels,[i,j])), 'snr'))');
            best_snr2 = max(compute_snr(selected_data_gene_cluster_mean(:,ismember(cell_labels,i)), selected_data_gene_cluster_mean(:,ismember(cell_labels,j))));
            best_snr3 = compute_snr_between_within(selected_data_gene_cluster_mean(:,ismember(cell_labels,i)),selected_data_gene_cluster_mean(:,ismember(cell_labels,j)));
            tmp = [tmp; [i,j,best_snr, best_snr2, best_snr3]];
        end
    end
    if isempty(tmp) || min(tmp(:,3))>snr_merge_threshold
        break;
    else
        [~,merge_entry_ind] = min(tmp(:,3));
        cell_labels(cell_labels==tmp(merge_entry_ind,2)) = tmp(merge_entry_ind,1);
        cell_labels = standardize_idx(cell_labels);
    end
end
fprintf(1, 'Cell-cell graph contains %d cell types after merging\n', length(unique(cell_labels)));
cell_cluster_best_snr = tmp;


% reorder the cell types and pathways for visualization
if max(cell_labels)==1
    return
end
pathway_in_cell_detection_rate = grpstats(data(gene_labels(keep_ind)~=0,:), gene_labels(keep_ind(gene_labels(keep_ind)~=0)));
pathway_in_celltype_detection_rate = grpstats(pathway_in_cell_detection_rate',cell_labels')';
pathway_in_celltype_detection_rate_normalize = pathway_in_celltype_detection_rate./max(pathway_in_celltype_detection_rate,[],2);
[U,S,V] = svds(pathway_in_celltype_detection_rate_normalize,1);
[~,I_celltype] = sort(V);
cell_labels_new = zeros(size(cell_labels));
for i=1:length(I_celltype)
    cell_labels_new(cell_labels==I_celltype(i)) = i;
end
cell_labels = cell_labels_new;

useful = zeros(size(pathway_in_cell_detection_rate,1),1);
scores = zeros(size(pathway_in_cell_detection_rate,1),max(cell_labels));
for k=1:max(cell_labels)
    for l=k+1:max(cell_labels)
        tmp = get_correlations(pathway_in_cell_detection_rate(:,ismember(cell_labels,[k,l])),cell_labels(ismember(cell_labels,[k,l])),'snr');
        useful = useful + double(abs(tmp)==max(abs(tmp)) | abs(tmp) >=snr_merge_threshold);
        scores(:,k) = scores(:,k) - tmp;
        scores(:,l) = scores(:,l) + tmp;
    end
end
[max_score,signature_for_which_celltype] = max(scores,[],2);
pathway_order = [];
for k=1:max(cell_labels)
    tmp_pathway = find(signature_for_which_celltype==k);
    [~,I] = sort(max_score(tmp_pathway),'descend');
    pathway_order = [pathway_order;tmp_pathway(I)];
end    
pathway_order(ismember(pathway_order, find(useful==0)))=[];
signature_for_which_celltype = signature_for_which_celltype(pathway_order);
gene_labels_new = zeros(size(gene_labels));
for i=1:length(pathway_order)
    gene_labels_new(gene_labels==pathway_order(i)) = i;
end
gene_labels = gene_labels_new;
fprintf(1, 'Number of useful pathways is %d\n', max(gene_labels));
