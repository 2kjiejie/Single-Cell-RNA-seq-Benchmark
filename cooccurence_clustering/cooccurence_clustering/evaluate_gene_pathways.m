function [snr_score, mena_diff_score, mean_ratio_score] = evaluate_gene_pathways(data1,data2)
% each colume is a cell, each row is a pathway

mena_diff_score = abs(mean(data1,2)-mean(data2,2));
mean_ratio_score = exp(abs(log( (mean(data1,2)+(1e-10)) ./ (mean(data2,2)+(1e-10)) )));       
snr_score = abs(get_correlations([data1, data2], [zeros(1,size(data1,2)),ones(1,size(data2,2))], 'snr'));



