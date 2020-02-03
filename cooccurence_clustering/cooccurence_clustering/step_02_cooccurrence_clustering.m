%% addpath to all tools
addpath(genpath('../cooccurence_clustering/'))


%% initiate one instance of the "cooccurrence_clustering_analysis" class
cooc = cooccurrence_clustering_analysis;


%% read data prepawed in matlab file
cooc = cooc.Read10X('C:\Users\jzhou417\Desktop\pbmc3k\cooc1');
%cooc = cooc.ReadMatlab('GSM3022245_gene_count.dat.csv');


%% filter the data by removing genes and cells (same as Seurat tutorial on this data)
cooc.initial_filtering_min_num_cells = 10;
cooc.initial_filtering_min_num_genes = 0;
cooc.initial_filtering_max_num_genes = Inf;
cooc.initial_filtering_max_percent_mito = 1;
cooc = cooc.initial_filtering_of_data(1);


%% binarize data
cooc.binarization_threshold = 0;
cooc.binary_data = full(double(cooc.data>cooc.binarization_threshold));



%% cooccurrence clustering
cooc.cooccurrence_min_expressed_cells = 10;             % genes will only be considered if detected in >= minimum number of cells, and undetected in >= minimum number of cells
cooc.cooccurrence_min_pathway_size = 20;             % only considered gene clusters of size >= this threshold
cooc.cooccurrence_min_population_size = 10;          % only considered gene clusters of size >= this threshold
cooc.cooccurrence_snr_merge_threshold = 1.5;            % threshold for merging Louvain communities, based on snr of average detection of each gene cluster
cooc.cooccurrence_mean_diff_merge_threshold = 0.5;      % threshold for merging Louvain communities
cooc.cooccurrence_mean_ratio_merge_threshold = 2;       % threshold for merging Louvain communities


cooc = cooc.iterative_cooccurance_clustering;


%% re-order pathways for better visualization

pathway_in_cell_detection_rate = cooc.pathway_in_cell_detection_rate;
scores = zeros(size(pathway_in_cell_detection_rate,1),max(cooc.cell_labels));
for k=1:max(cooc.cell_labels)
    for l=k+1:max(cooc.cell_labels)
        tmp = get_correlations(pathway_in_cell_detection_rate(:,ismember(cooc.cell_labels,[k,l])),cooc.cell_labels(ismember(cooc.cell_labels,[k,l])),'snr');
        scores(:,k) = scores(:,k) - tmp;
        scores(:,l) = scores(:,l) + tmp;
    end
end
[max_score,signature_for_which_celltype] = max(scores,[],2);
pathway_order = [];
for k=1:max(cooc.cell_labels)
    tmp_pathway = find(signature_for_which_celltype==k);
    [~,I] = sort(max_score(tmp_pathway),'descend');
    pathway_order = [pathway_order;tmp_pathway(I)];
end    
cooc.pathway_in_cell_detection_rate = cooc.pathway_in_cell_detection_rate(pathway_order,:);
cooc.gene_pathway_indicator = cooc.gene_pathway_indicator(:,pathway_order);
                    
save([mfilename,'_result.mat']);



%% Visualize the percentages of detecion for pathways in individual cells

[~,I] = sort(cooc.cell_labels);
h = figure(30); set(h, 'position', [100 100 1200 600]);
subplot(6,1,1:4); imagesc(cooc.pathway_in_cell_detection_rate(:, I)./max(cooc.pathway_in_cell_detection_rate(:, I),[],2)); ylabel('pathways (gene clusters)'); xlabel('cells'); title('Percentage of detection')
h = subplot(4,1,4); imagesc(cooc.cell_labels(I)); h.YTick=[];
for k=1:max(cooc.cell_labels)
    text(mean(find(cooc.cell_labels(I)==k)), 1, num2str(unique(cooc.cell_labels_history(cooc.cell_labels==k,end))));
end
title('co-occurrence cell clusters'); xlabel('cells')


h = figure(40); set(h, 'position', [100 100 1200 600]);
subplot(6,1,1:4); imagesc(cooc.pathway_in_cell_detection_rate(:, I)./max(cooc.pathway_in_cell_detection_rate(:, I),[],2)); ylabel('pathways (gene clusters)'); xlabel('cells'); title('Percentage of detection')
h = subplot(4,1,4); imagesc(cooc.cell_labels(I)); h.YTick=[];
for k=1:max(cooc.cell_labels)
    text(mean(find(cooc.cell_labels(I)==k)), 1, num2str(k));
end
title('co-occurrence cell clusters (renamed)'); xlabel('cells')



%%

h = figure(5); set(h, 'units','normalized','outerposition',[0 0 1 1])
for i=1:size(cooc.pathway_in_cell_detection_rate,1)
    tmp = cell(0);
    for k=1:max(cooc.cell_labels)
        tmp{k} = cooc.pathway_in_cell_detection_rate(i,cooc.cell_labels==k);
    end
    subplot(5, ceil(size(cooc.pathway_in_cell_detection_rate,1)/5),i); violin(tmp,'mc',[],'medc',[]);
    title(num2str(i,'Pathway - %d'))
    xlabel('cell clusters');
    ylabel('% detection');
    drawnow
end



%% export results to files

pathway_names = strcat('pathway-',regexp(num2str(1:size(cooc.gene_pathway_indicator,2),'%05d,'),',','split')'); pathway_names(end)=[];
%write_to_txt_v2([mfilename,'_result_gene_pathway_indicator.csv'],[{' '}, pathway_names'], cooc.gene_names, cooc.gene_pathway_indicator,',');
%write_to_txt_v2([mfilename,'_result_pathway_in_cell_detection_rate.csv'],[{' '}, cooc.cell_names], pathway_names, cooc.pathway_in_cell_detection_rate,',');
write_to_txt_v2([mfilename,'_result_cell_clusters_pbmc3k_stan.csv'],[{' '}, cooc.cell_names], {'cluster idx'}, cooc.cell_labels,',');

save('results_pbmc3k_stan.mat')


