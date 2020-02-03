function [idx] = clustering_graph_by_modularity(graph_file, resolution)

if ~exist('resolution')
    resolution = 1;
end

fprintf('Running ModularityOptimizer for clustering ...');
tmp = regexp(path,':','split');
% if sum(path=='\')>sum(path=='/')
%     seurat_java_path = tmp{union(isInListEnd(tmp,'eurat\java'), isInListEnd(tmp,'modularity_java'))};
% else
%     seurat_java_path = tmp{union(isInListEnd(tmp,'eurat/java'), isInListEnd(tmp,'modularity_java'))};
% end
%input the jave path as the optimizer
seurat_java_path = 'C:\Users\jzhou417\Desktop\cooccurence_clustering\cooccurence_clustering\Seurat\java';
tic
system(['java -jar ', fullfile(seurat_java_path,'ModularityOptimizer.jar'), ...
        ' ', fullfile(cd, graph_file), ...
        ' ', fullfile(cd, 'output_file.txt'), ...
        ' 1', ...                           % Modularity function to use in clustering (1 = standard; 2 = alternative).
        ' ', num2str(resolution), ...       % Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
        ' 2', ...                           % Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
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
idx = str2double(cluster_names)+1;
delete(fullfile(cd, 'output_file.txt'));
delete(fullfile(cd, graph_file));


