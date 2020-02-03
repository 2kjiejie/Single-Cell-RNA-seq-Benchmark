function [cell_labels_history, gene_labels_history,cell_labels, standardized_cell_labels, updated_standardized_cell_labels, gene_pathway_indicator] = iterative_cooccurance_clustering(binary_data, min_expressed_cells, min_cluster_genes, snr_merge_threshold)

rng(500);

cell_labels_history = [zeros(size(binary_data,2),1)];
gene_labels_history = [zeros(size(binary_data,1),1)];
cell_labels = [zeros(size(binary_data,2),1)];
cell_clusters_to_process = 0;
while ~isempty(cell_clusters_to_process)
    fprintf(1,'Remaining clusters to partition %d\n', length(cell_clusters_to_process));
    cells_considered = find(cell_labels == cell_clusters_to_process(1));
    [cell_labels_tmp, gene_labels_tmp, cell_cluster_best_snr] = cooccurance_clustering(binary_data(:,cells_considered), min_expressed_cells, min_cluster_genes, snr_merge_threshold);
    if length(unique(cell_labels_tmp))==1
        cell_clusters_to_process(1)=[];
    else
        cell_cluster_best_snr
        
        h = figure(1); set(h,'Position',[100,50,1000,700]);
        selected_data = binary_data(gene_labels_tmp~=0,cells_considered);
        selected_gene_labels = gene_labels_tmp(gene_labels_tmp~=0);
        selected_cell_labels = cell_labels_tmp;
        [~,perm_r] = sort(selected_gene_labels);
        [~,perm_c] = sort(selected_cell_labels);
        subplot(2,4,[1 2 3]);imagesc(selected_data(perm_r, perm_c))
        subplot(2,4,[4]);imagesc([selected_gene_labels(perm_r)])
        subplot(6,4,[13:15]);imagesc([selected_cell_labels(perm_c)'])
        subplot(6,4,[17:19 21:23]);imagesc(grpstats(selected_data(:,perm_c), selected_gene_labels));
        drawnow
        snapnow

        cell_clusters_to_process(1)=[];
        cell_clusters_to_process = [cell_clusters_to_process; unique(max(cell_labels)+cell_labels_tmp)];
        cell_labels(cells_considered) = max(cell_labels)+cell_labels_tmp;
        cell_labels_history = [cell_labels_history, cell_labels];
        gene_labels_history = [gene_labels_history, gene_labels_tmp];
    end
    fprintf(1,'\n');
end



% standardize the cell labels and identify useful pathways
standardized_cell_labels = standardize_idx(cell_labels)';

pathway_in_cell_detection_rate = [];
gene_pathway_indicator = [];
for i=1:size(gene_labels_history,2)
    for j=1:max(gene_labels_history(:,i))
        gene_pathway_indicator = [gene_pathway_indicator, gene_labels_history(:,i)==j];
        pathway_in_cell_detection_rate = [pathway_in_cell_detection_rate; mean(binary_data(gene_labels_history(:,i)==j,:)~=0,1)];
    end
end

useful = zeros(size(pathway_in_cell_detection_rate,1),max(standardized_cell_labels));
for k=1:max(standardized_cell_labels)
    for l=k+1:max(standardized_cell_labels)
        tmp = abs(get_correlations(pathway_in_cell_detection_rate(:,ismember(standardized_cell_labels,[k,l])),standardized_cell_labels(ismember(standardized_cell_labels,[k,l])),'snr'));
        % sum(tmp >=snr_merge_threshold)
        useful(:,k) = useful(:,k) + double(tmp==max(tmp) | tmp >=snr_merge_threshold);
        useful(:,l) = useful(:,l) + double(tmp==max(tmp) | tmp >=snr_merge_threshold);
    end
end


% use gaussian to model each cell label
mu=[];
sigma = [];
for i=1:max(standardized_cell_labels)
    mu = [mu, mean(pathway_in_cell_detection_rate(:,standardized_cell_labels==i),2)];
    sigma = [sigma, std(pathway_in_cell_detection_rate(:,standardized_cell_labels==i),[],2)];
end
sigma(sigma==0)=min(sigma(sigma~=0));

% mapping the cells
prob = ones(max(standardized_cell_labels), size(pathway_in_cell_detection_rate,2));
for i=1:max(standardized_cell_labels)
    for j=1:size(gene_pathway_indicator,2)
        if useful(j,i)~=0
            prob(i,:) = prob(i,:).*normpdf(pathway_in_cell_detection_rate(j,:),mu(j,i),sigma(j,i));
        end
    end
end
prob = prob./sum(prob,1);

% update labels
[Y,updated_standardized_cell_labels] = max(prob,[],1);
updated_standardized_cell_labels(Y<0.8)=0;
updated_standardized_cell_labels(updated_standardized_cell_labels~=standardized_cell_labels)=0;
