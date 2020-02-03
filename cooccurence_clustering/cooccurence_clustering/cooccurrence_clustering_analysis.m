classdef cooccurrence_clustering_analysis
    properties
        raw_data                % number of parameters
        raw_gene_names          % number of dynamic variables
        raw_cell_names          % number of dynamic variables

        genes_to_keep
        cells_to_keep
        
        data
        binary_data
        cell_names
        gene_names
        
        high_variance_genes
        residue_data

        PCA_data_used
        PCA_gene_used
        PCA_cell_embeddings
        PCA_gene_loadings
        PCA_gene_loadings_full
        PCA_variances
        jackstraw_random_loadings
        jackstraw_real_loadings
        
        meta_data
        meta_feature_names      
        
        pathway_in_cell_detection_rate
        gene_pathway_indicator
        gene_labels_history
        cell_labels_history
        vis_ordering_cell_labels
        vis_ordering_gene_labels
        cell_labels
        
        cell_labels_seurat

        cooccurrence_min_expressed_cells = 10;        % genes will only be considered if detected in >= minimum number of cells, and undetected in >= minimum number of cells
        cooccurrence_min_pathway_size = 10;           % only considered gene clusters of size >= this threshold
        cooccurrence_min_population_size = 10;        % cell clusters cannot be smaller than this threshold, otherwise will be merged with other clusters
        cooccurrence_snr_merge_threshold = 1.5;       % threshold for merging Louvain communities, based on snr of average detection of each gene cluster
        cooccurrence_mean_diff_merge_threshold = 0.5;       % threshold for merging Louvain communities
        cooccurrence_mean_ratio_merge_threshold = 2;       % threshold for merging Louvain communities

        nUMI_scaling_factor = 10000;                % scale_cells_by_nUMI nUMI factor
        normalization_method = 'none';              % 'none' or 'log1p'
        binarization_threshold = 0;                 % threshold for turning data into binary
        
        nBins_for_HVG = 20;                         % number of bins for computing high variance genes
        expmean_threshold_low = 0.0125;             % expmean threshold for selecting HVG
        expmean_threshold_high = 3;                 % expmean threshold for selecting HVG
        dispersion_threshold_low = 0.5;             % dispersion threshold for selecting HVG
        
        initial_filtering_min_num_cells = 0;        % initial filtering of the raw data
        initial_filtering_min_num_genes = 0;
        initial_filtering_max_num_genes = Inf;
        initial_filtering_max_percent_mito = 1;
        
        num_PC = 20;
    end
    
    methods
        %%
        function obj = cooccurrence_clustering_analysis  % initiate class
            obj.raw_data = [];
            obj.raw_gene_names = [];
            obj.raw_cell_names = [];
        end       
        
        %%
        function obj = Read10X(obj, folder_name)
            
            % cell names
            fprintf('Loading cell names ... ');
            tic
            fid = fopen(fullfile(folder_name,'barcodes.tsv'));
            d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
            d(d==char(13))=[];
            cell_names = regexp(d,char(10),'split');
            cell_names(end)=[];
            obj.raw_cell_names = cell_names;
            toc
            % gene names
            fprintf('Loading gene names ... ');
            tic
            fid = fopen(fullfile(folder_name,'genes.tsv'));
            d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
            d(d==char(13))=[];
            d(d==char(9)) = char(10);
            gene_names = regexp(d,char(10),'split');
            %gene_names = gene_names(2:2:end)';
            
            
            obj.raw_gene_names = gene_names;
            toc
            % raw_data
            fprintf('Loading raw data ... ');
            tic
            fid = fopen(fullfile(folder_name,'matrix.mtx'));
            d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
            d(1:max(find(d==char(10),2)))=[];
            d(d==char(10))=char(9);
            s = str2num(d);
            num_genes = s(1);
            num_cells = s(2);
            num_entries = s(3);
            raw = reshape(s(4:end),3,length(s)/3-1);
            toc
            % reformat
            fprintf('Converting data of %d genes * %d cells into matrix ... ',num_genes, num_cells);
            tic
            obj.raw_data = sparse(raw(1,:),raw(2,:), raw(3,:));
            if size(obj.raw_data,1)~=num_genes || size(obj.raw_data,2)~=num_cells
                obj.raw_data(num_genes,num_cells)=0;
            end
            toc
        end        
        
        %%
        function obj = ReadMatlab(obj, filename)
            %load(filename,'raw_data','raw_gene_names','raw_cell_names');
            T = readtable(filename,'ReadRowNames',true);
            obj.raw_data = table2array(T);
            obj.raw_gene_names = T.Properties.RowNames;
            obj.raw_cell_names = T.Properties.VariableNames;
        end
        
        %%
        function obj = initial_filtering_of_data(obj, isdisplay)
            if ~exist('isdisplay')
                isdisplay = 0;
            end
            
            nGene = full(sum(obj.raw_data~=0,1));
            nUMI = full(sum(obj.raw_data,1));
            mito_genes = obj.raw_gene_names(isInListFront(obj.raw_gene_names,'MT-'));
            percent_mito = full(sum(obj.raw_data(ismember(obj.raw_gene_names, mito_genes),:),1)./sum(obj.raw_data,1));
            
            if isdisplay==1
                figure(1); subplot(2,3,1); violin(nGene'); title('nGene')
                figure(1); subplot(2,3,2); violin(nUMI'); title('nUMI')
                figure(1); subplot(2,3,3); violin(percent_mito'); title('percent mito')
                figure(1); subplot(2,2,3); plot(nUMI, percent_mito, '.'); xlabel('nUMI'); ylabel('percent mito')
                figure(1); subplot(2,2,4); plot(nUMI, nGene, '.'); xlabel('nUMI'); ylabel('nGene')
            end
            
            obj.genes_to_keep = find(sum(obj.raw_data>0,2) >= obj.initial_filtering_min_num_cells);    
            obj.cells_to_keep = find(nGene>=obj.initial_filtering_min_num_genes & nGene<=obj.initial_filtering_max_num_genes & percent_mito<=obj.initial_filtering_max_percent_mito); 
            
            obj.gene_names = obj.raw_gene_names(obj.genes_to_keep);
            obj.cell_names = obj.raw_cell_names(obj.cells_to_keep);
            obj.data = obj.raw_data(obj.genes_to_keep, obj.cells_to_keep);
            fprintf('Data after initial filtering %d genes * %d cells. \n', length(obj.gene_names), length(obj.cell_names));
        end
        
        %%
        function obj = scale_cells_by_nUMI(obj)
            try
                obj.data = obj.data./sum(obj.data,1)*obj.nUMI_scaling_factor;  
            catch
                obj.data = bsxfun(@rdivide,obj.data,sum(obj.data,1))*obj.nUMI_scaling_factor; 
            end
        end       
        
        %%
        function obj = normalize_data(obj)
            if isequal(obj.normalization_method,'log1p')
                obj.data = log(obj.data+1); 
                if sum(obj.data(:)~=0)<=0.2*prod(size(obj.data)) % data is still sparse
                    obj.data = sparse(obj.data); 
                end
            elseif isequal(obj.normalization_method,'log1p')
                % do nothing
            else
                fprintf('obj.normalization_method should be either "none" or "log1p"\n'); % 'none' or 'log1p'
            end
        end
        
        %% select high variant genes
        function obj = select_high_variance_genes(obj, isdisplay)
            if ~exist('isdisplay')
                isdisplay = 0;
            end
            expmean = log(mean(exp(obj.data)-1,2)+1);
            logvmr = log(var(exp(obj.data)-1,[],2)./mean(exp(obj.data)-1,2));
            logvmr(isnan(logvmr))=0;
            logvmr(isinf(logvmr))=0;
            [~,bin_centers] = hist(expmean,obj.nBins_for_HVG);  % parameter 20 bins
            try 
                [~,bin_assignment] = min(abs(expmean-bin_centers),[],2);
            catch
                [~,bin_assignment] = min(abs(bsxfun(@minus, expmean,bin_centers)),[],2);
            end
            
            scaled_logvmr = zeros(size(logvmr));
            for i=1:length(bin_centers)
                mu_tmp = mean(logvmr(bin_assignment==i));
                std_tmp = std(logvmr(bin_assignment==i));
                scaled_logvmr(bin_assignment==i) = (logvmr(bin_assignment==i)-mu_tmp)/std_tmp;
            end
            scaled_logvmr(isnan(scaled_logvmr) | isinf(scaled_logvmr))=0;
            obj.high_variance_genes = obj.gene_names(expmean>obj.expmean_threshold_low & expmean<obj.expmean_threshold_high & scaled_logvmr>=obj.dispersion_threshold_low);
            if isdisplay==1
                figure(1); subplot(1,1,1)
                plot(expmean, scaled_logvmr,'go',expmean(ismember(obj.gene_names,obj.high_variance_genes)), scaled_logvmr(ismember(obj.gene_names,obj.high_variance_genes)),'bo')
            end
            fprintf('Number of high variance genes is %d \n', length(obj.high_variance_genes));
        end
        
        %%
        function obj = add_metadata(obj,new_meta_data, new_feature)
            obj.meta_data = [obj.meta_data; new_meta_data];
            obj.meta_data = full(obj.meta_data);
            if ~iscell(new_feature)
                new_feature = {new_feature};
            end
            obj.meta_feature_names = [obj.meta_feature_names; new_feature];
        end
        
        %%
        function obj = regress_out_unwanted_variation(obj, features_to_regress_out)
            if ~exist('features_to_regress_out')
                features_to_regress_out = obj.meta_feature_names;
            end
            
            meta_data = [ones(1,size(obj.data,2)); obj.meta_data(ismember(obj.meta_feature_names, features_to_regress_out),:)];
            res = (obj.data' - meta_data'*(pinv(meta_data')*obj.data'))';
            try
                scaled_res = (res-mean(res,2))./std(res,[],2);
            catch
                scaled_res = bsxfun(@rdivide, (bsxfun(@minus,res,mean(res,2))), std(res,[],2));
            end
            scaled_res(scaled_res>10)=10; 
            scaled_res(isnan(scaled_res))=0;
            obj.residue_data = scaled_res;
        end    
        
        %%
        function obj = PCA(obj, data_used, genes_used)
            if ~exist('data_used')
                if ~isempty(obj.residue_data)
                    data_used = obj.residue_data;
                else
                    data_used = obj.data;
                end
            end
            if ~exist('genes_used')
                if ~isempty(obj.high_variance_genes)
                    genes_used = obj.high_variance_genes;
                else
                    genes_used = obj.gene_names;
                end
            end
            obj.PCA_data_used = data_used; 
            obj.PCA_gene_used = genes_used; 
            [U,S,V] = svds(data_used(ismember(obj.gene_names,genes_used),:),obj.num_PC);
            obj.PCA_cell_embeddings = V(:,:)*S;
            obj.PCA_gene_loadings = U(:,:);           
            obj.PCA_gene_loadings_full = data_used*(V);
            obj.PCA_variances = diag(S);
        end
    
        %%
        function [num_significant_PCs, max_random_variance] = PCA_significant_PC_by_permutation(obj, perm_iter, isdisplay)
            if ~exist('perm_iter')
                perm_iter = 10; 
            end
            if ~exist('isdisplay')
                isdisplay = 0;
            end
            
            data_used = obj.PCA_data_used(ismember(obj.gene_names, obj.PCA_gene_used),:);
            highest_random_variance = zeros(1,perm_iter);
            fprintf('Random permutation %d iter for computing signficant PC ... %5d ', perm_iter, 0);
            tic
            for i=1:perm_iter
                perm_data = zeros(size(data_used));
                [~,I] = sort(rand(size(data_used)),2);
                perm_data = data_used((I-1)*size(data_used,1)+repmat((1:size(data_used,1))',1,size(data_used,2)));
                [U,S,V] = svds(perm_data,obj.num_PC);
                highest_random_variance(i) = max(diag(S));
                fprintf('\b\b\b\b\b\b%5d ', i);
            end
            toc
            max_random_variance = max(highest_random_variance);
            num_significant_PCs = sum(obj.PCA_variances>=max_random_variance);
            if isdisplay==1
                plot(obj.PCA_variances,'o')
                line(xlim,max(highest_random_variance)*[1 1])
                title(num2str(num_significant_PCs, 'Significant PCs = %d')); 
                xlabel('PC'); ylabel('standard deviation')
            end
        end
        
        %%
        function [num_significant_PCs, pvalue] = PCA_significant_PC_by_jackstraw(obj, perm_iter, proportion_to_permute, isdisplay)
            if ~exist('perm_iter')
                perm_iter = 10; 
            end
            if ~exist('proportion_to_permute')
                proportion_to_permute = 0.01; 
            end
            if ~exist('isdisplay')
                isdisplay = 0;
            end
            
            data_used = obj.PCA_data_used(ismember(obj.gene_names, obj.PCA_gene_used),:);
            fprintf('Jackstraw %d iter for computing signficant PC ... %5d ', perm_iter, 0);
            tic
            jackstraw_real_loadings = [];
            jackstraw_random_loadings = [];
            for i=1:perm_iter
                perm_data = data_used;
                ind_to_permute = sort(randsample(1:size(perm_data,1),round(size(perm_data,1)*proportion_to_permute)));
                [~,I] = sort(rand(length(ind_to_permute), size(perm_data,2)),2);
                perm_data(ind_to_permute,:) = perm_data((I-1)*size(perm_data,1)+repmat(ind_to_permute(:),1,size(perm_data,2)) );
                [U,S,V] = svds(perm_data,obj.num_PC);
                jackstraw_real_loadings = [jackstraw_real_loadings; obj.PCA_gene_loadings(ind_to_permute,:)];
                jackstraw_random_loadings = [jackstraw_random_loadings; U(ind_to_permute,:)];
                fprintf('\b\b\b\b\b\b%5d ', i);
            end
            toc
            [~,pvalue] = ttest(abs(jackstraw_real_loadings),abs(jackstraw_random_loadings),'tail','right');
            num_significant_PCs = min(find(pvalue>0.01/obj.num_PC))-1;  % PCs are significant up until the first one that isn't
            if isempty(num_significant_PCs)
                num_significant_PCs = obj.num_PC;
            end
            obj.jackstraw_real_loadings = jackstraw_real_loadings;
            obj.jackstraw_random_loadings = jackstraw_random_loadings;
            if isdisplay==1
                h = figure;
                set(h, 'units','normalized','outerposition',[0 0 1 1])
                for i=1:obj.num_PC
                    subplot(4,ceil(obj.num_PC/4),i);
                    plot(sort(abs(jackstraw_random_loadings(:,i))),sort(abs(jackstraw_real_loadings(:,i))),'o');
                    axis([0,1,0,1]*max(abs(axis)))
                    line(xlim,ylim)
                    title(num2str([i, pvalue(i)],'PC %d - pvalue: %f'));
                end
            end
        end
        
        %%
        function obj = iterative_cooccurance_clustering(obj)
            
            rng(500);
            cell_labels_history = [zeros(size(obj.binary_data,2),1)];
            gene_labels_history = [zeros(size(obj.binary_data,1),1)];
            cell_labels = [zeros(size(obj.binary_data,2),1)];
            cell_clusters_to_process = 0;
            
            while ~isempty(cell_clusters_to_process)
                fprintf(1,'Remaining clusters to partition %d\n', length(cell_clusters_to_process));
                fprintf(1,'Processing cluster %d now ...\n', cell_clusters_to_process(1));
                cells_considered = find(cell_labels == cell_clusters_to_process(1));
                [cell_labels_tmp, gene_labels_tmp, cell_cluster_best_snr, signature_for_which_celltype] = cooccurance_clustering(obj.binary_data(:,cells_considered), obj.cooccurrence_min_expressed_cells, obj.cooccurrence_min_pathway_size, obj.cooccurrence_min_population_size, obj.cooccurrence_snr_merge_threshold, obj.cooccurrence_mean_diff_merge_threshold, obj.cooccurrence_mean_ratio_merge_threshold);

                if length(unique(cell_labels_tmp))==1
                    cell_clusters_to_process(1)=[];
                else
                    % cell_cluster_best_snr
                    
                    
                    h = figure(1); set(h,'Position',[100,50,1000,700]);
                    selected_data = obj.binary_data(gene_labels_tmp~=0,cells_considered);
                    selected_gene_labels = gene_labels_tmp(gene_labels_tmp~=0);
                    selected_cell_labels = cell_labels_tmp;
                    [~,perm_r] = sort(selected_gene_labels+rand(size(selected_gene_labels))/10);
                    [~,perm_c] = sort(selected_cell_labels+rand(size(selected_cell_labels))/10);
                    subplot(5,7,setdiff(1:3*7, (1:3)*7));imagesc(selected_data(perm_r, perm_c)); title(['cluster ', num2str(unique(cell_labels_history(cells_considered,end)))],'fontsize',15); ylabel('selected genes','fontsize',15)
                    h = subplot(5,9,[(1:3)*9]);imagesc([selected_gene_labels(perm_r)]); h.XTick=[]; title('gene pathways','fontsize',15)
                    for k=1:max(selected_gene_labels)
                        text(1,mean(find(selected_gene_labels(perm_r)==k)), char('a'-1+k),'fontsize',15);
                    end
                    h = subplot(13,7,[64:69]);imagesc([selected_cell_labels(perm_c)']); h.YTick=[]; h.XTick=[]; title('cell subclusters','fontsize',15)
                    for k=1:max(selected_cell_labels)
                        text(mean(find(selected_cell_labels(perm_c)==k)), 1, num2str(k+max(cell_labels)),'fontsize',15);
                    end
                    try
                        h = subplot(5,7,[29:34]);imagesc(grpstats(selected_data(:,perm_c), selected_gene_labels)./max(grpstats(selected_data(:,perm_c), selected_gene_labels),[],2));
                    catch
                        h = subplot(5,7,[29:34]);imagesc(bsxfun(@rdivide, grpstats(selected_data(:,perm_c), selected_gene_labels), max(grpstats(selected_data(:,perm_c), selected_gene_labels),[],2)));
                    end
                    h.YTick=1:max(selected_gene_labels);
                    h.YTickLabel = mat2cell(char((1:max(selected_gene_labels))+'a'-1)',ones(max(selected_gene_labels),1),1);
                    xlabel('cells','fontsize',15);
                    ylabel('gene pathways','fontsize',15);
                    title('Pathway activities','fontsize',15)
                    drawnow
                    snapnow

                    % figure out how to order the cell types and pathways in the final combined visualization
                    tmp_vis_ordering_cell_labels = [1:max(cell_labels_tmp);1:max(cell_labels_tmp)];
                    if isempty(setdiff(1:max(cell_labels_tmp), signature_for_which_celltype))
                        tmp_vis_ordering_gene_labels = [1:max(gene_labels_tmp);(1:max(gene_labels_tmp));signature_for_which_celltype(:)'];
                    else
                        tmp_vis_ordering_gene_labels = [1:max(gene_labels_tmp);(1:max(gene_labels_tmp));signature_for_which_celltype(:)'];
                        tmp2 = [];
                        for i=1:max(cell_labels_tmp)
                            ind = find(tmp_vis_ordering_gene_labels(3,:)==i);
                            if isempty(ind)
                                tmp2 = [tmp2,[0;NaN;i]];
                            else
                                tmp2 = [tmp2,tmp_vis_ordering_gene_labels(:,ind)];
                            end
                        end
                        tmp2(2,isnan(tmp2(2,:))) = interp1([0, setdiff(1:length(tmp2(2,:)),find(isnan(tmp2(2,:)))), length(tmp2(2,:))+1], [0, tmp2(2,~isnan(tmp2(2,:))), length(tmp2(2,:))+1], find(isnan(tmp2(2,:))), 'linear', 'extrap');
                        tmp2(2,:) = tmp2(2,:) - min(tmp2(2,:)) + 1;
                        tmp_vis_ordering_gene_labels = tmp2;
                    end
                    
                    if isequal(cells_considered, (1:length(cell_labels))')
                        vis_ordering_cell_labels = tmp_vis_ordering_cell_labels;
                        vis_ordering_gene_labels = tmp_vis_ordering_gene_labels;
                    else
                        mother_celltype_vis_position = vis_ordering_cell_labels(2,ismember(vis_ordering_cell_labels(1,:), cell_clusters_to_process(1)));
                        next_celltype_vis_position = min(vis_ordering_cell_labels(2,vis_ordering_cell_labels(2,:)>mother_celltype_vis_position));
                        if isempty(next_celltype_vis_position)
                            next_celltype_vis_position = mother_celltype_vis_position + max(cell_labels_tmp) + 1;
                        end
                        celltype_vis_positions = interp1([0,max(cell_labels_tmp)+1], [mother_celltype_vis_position,next_celltype_vis_position], 1:max(cell_labels_tmp));
                        vis_ordering_cell_labels = [vis_ordering_cell_labels,...
                                                    [(1:max(cell_labels_tmp)) + max(cell_labels); ...
                                                     celltype_vis_positions] ...
                                                   ];
                        
                        mother_celltype_last_pathway_vis_position = max(vis_ordering_gene_labels(2,ismember(vis_ordering_gene_labels(3,:), cell_clusters_to_process(1))));
                        next_celltype_first_pathway_vis_position = min(vis_ordering_cell_labels(2,vis_ordering_cell_labels(2,:)>mother_celltype_last_pathway_vis_position));
                        if isempty(next_celltype_first_pathway_vis_position)
                            next_celltype_first_pathway_vis_position = mother_celltype_last_pathway_vis_position + max(signature_for_which_celltype) + 1;
                        end
                        pathway_vis_positions = interp1([0,length(tmp_vis_ordering_gene_labels(2,:))+1], [mother_celltype_last_pathway_vis_position,next_celltype_first_pathway_vis_position],tmp_vis_ordering_gene_labels(2,:));
                        vis_ordering_gene_labels = [vis_ordering_gene_labels, ...
                                                    [tmp_vis_ordering_gene_labels(1,:) + (tmp_vis_ordering_gene_labels(1,:)~=0).*max(max(gene_labels_history)); ...
                                                     pathway_vis_positions; ...
                                                     tmp_vis_ordering_gene_labels(3,:) + max(cell_labels)] ...
                                                    ];
                    end
                    
                    cell_clusters_to_process(1)=[];
                    cell_clusters_to_process = [cell_clusters_to_process; unique(max(cell_labels)+cell_labels_tmp)];
                    cell_labels(cells_considered) = max(cell_labels)+cell_labels_tmp;
                    cell_labels_history = [cell_labels_history, cell_labels];
                    gene_labels_history = [gene_labels_history, gene_labels_tmp + (gene_labels_tmp~=0)*max(max(gene_labels_history))];
                end
                fprintf(1,'\n');
            end


            % standardize the cell labels and identify useful pathways
            if ~exist('vis_ordering_cell_labels') || isempty(vis_ordering_cell_labels)
                standardized_cell_labels = cell_labels;
                pathway_in_cell_detection_rate = [];
                gene_pathway_indicator = [];
            else
                standardized_cell_labels = zeros(size(cell_labels));
                [~,I] = sort(vis_ordering_cell_labels(2,:));
                counter = 1;
                for i=1:length(I)
                    if sum(cell_labels==vis_ordering_cell_labels(1,I(i)))~=0
                        standardized_cell_labels(cell_labels==vis_ordering_cell_labels(1,I(i)))=counter;
                        counter = counter + 1;
                    end
                end

                [~,I] = sort(vis_ordering_gene_labels(2,:));
                gene_pathway_indicator=[];
                pathway_in_cell_detection_rate=[];
                for i=1:length(I)
                    if vis_ordering_gene_labels(1,I(i))==0
                        continue;
                    end
                    gene_pathway_indicator = [gene_pathway_indicator, sum(gene_labels_history==vis_ordering_gene_labels(1,I(i)),2)~=0];
                    pathway_in_cell_detection_rate = [pathway_in_cell_detection_rate; mean(obj.binary_data(gene_pathway_indicator(:,end)==1,:)~=0,1)];
                end
            end            
            
            obj.pathway_in_cell_detection_rate = pathway_in_cell_detection_rate;
            obj.gene_pathway_indicator = gene_pathway_indicator;
            obj.cell_labels_history = cell_labels_history;
            obj.gene_labels_history = gene_labels_history;
            obj.cell_labels = standardized_cell_labels(:)';

        end
        
        
        %%
        function obj = seurat_clustering(obj, data_used, SNN_k, SNN_prune, resolution)
            % data_used: each row is one cell, each column is one feature, i.e. cell_embeddings(:,1:10)
            % SNN_k: 30
            % SNN_prune: 1/15
            % resolution: 0.8
            
            
            % build a graph 
            [A, graph_file]=Seurat_SNN(data_used, SNN_k, SNN_prune);
            
            % modularity clustering
            fprintf('Running ModularityOptimizer for clustering ...');
            tmp = regexp(path,';','split');
            if sum(path=='\')>sum(path=='/')
                seurat_java_path = tmp{union(isInListEnd(tmp,'seurat\java'), isInListEnd(tmp,'modularity_java'))};
            else
                seurat_java_path = tmp{union(isInListEnd(tmp,'seurat/java'), isInListEnd(tmp,'modularity_java'))};
            end
            tic
            system(['java -jar ', fullfile(seurat_java_path,'ModularityOptimizer.jar'), ...
                    ' ', fullfile(cd, graph_file), ...
                    ' ', fullfile(cd, 'output_file.txt'), ...
                    ' 1', ...                           % Modularity function to use in clustering (1 = standard; 2 = alternative).
                    ' ', num2str(resolution), ...       % Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
                    ' 1', ...                           % Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
                    ' 100', ...                         % Number of random starts
                    ' 10', ...                          % Maximal number of iterations per random start.
                    ' 0', ...                           % Seed of the random number generator
                    ' 0', ...                           % Whether or not to print output to the console
                    ]);
            toc
            fid = fopen(fullfile(cd, 'output_file.txt'));
            d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
            d(d==char(13))=[];
            cluster_names = regexp(d,char(10),'split');
            cluster_names(end)=[];
            cell_idx = str2double(cluster_names)+1;
            delete(fullfile(cd, 'output_file.txt'));
            delete(fullfile(cd, graph_file));
            
            % merge singletons (not implemented yet)
            obj.cell_labels_seurat = cell_idx;
            
        end
        
    end
    
end


