function annotations = MSigDB_annotate_geneset(geneset, all_gene_names, MSigDB_categories_to_query)

if ~exist('MSigDB_categories_to_query')
    MSigDB_categories_to_query = {'h.all', 'c1.all', 'c2.cgp', 'c2.cp', 'c3.mir', 'c3.tft', 'c4.cgn', 'c4.cm', 'c5.bp', 'c5.cc', 'c5.mf', 'c6.all', 'c7.all'}
end    
  
load('MSigDB_processed_data.mat')

genesets_to_keep = ismember(MSigDB_geneset_category,MSigDB_categories_to_query);
MSigDB_geneset_category = MSigDB_geneset_category(genesets_to_keep);
MSigDB_geneset_membership = MSigDB_geneset_membership(:,genesets_to_keep);
MSigDB_geneset_name = MSigDB_geneset_name(genesets_to_keep);
MSigDB_geneset_ref = MSigDB_geneset_ref(genesets_to_keep);


[overlapping_genes,IA,~] = intersect(MSigDB_gene_symbol, all_gene_names);

pathhway_to_annotate = geneset;
total_size = length(overlapping_genes);
geneset_size = full(sum(MSigDB_geneset_membership(IA,:),1));
pathway_size = length(intersect(pathhway_to_annotate, overlapping_genes));
overlap_size = full(sum(MSigDB_geneset_membership(ismember(MSigDB_gene_symbol,pathhway_to_annotate),:),1));
pvalues = 1-hygecdf(overlap_size-1,total_size, geneset_size, pathway_size);
[Y,I] = sort(pvalues);

top_pathways_to_return = min(50, sum(pvalues<0.05));
annotations = [MSigDB_geneset_name(I(1:top_pathways_to_return)),MSigDB_geneset_category(I(1:top_pathways_to_return)), num2cell(pvalues(I(1:top_pathways_to_return)),1)'];
%[i,pathway_size, length(pathhway_to_annotate)]
