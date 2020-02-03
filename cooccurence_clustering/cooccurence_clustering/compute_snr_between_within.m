function r = compute_snr_between_within(data1, data2)

% each row is detection rate of one pathway
% each column is one cell

mean1 = mean(data1,2);
mean2 = mean(data2,2);
mean_all = (mean1+mean2)/2;

within1 = mean(sum((data1-repmat(mean1,1,size(data1,2))).^2,1));
within2 = mean(sum((data2-repmat(mean2,1,size(data2,2))).^2,1));

between = (mean(sum((data1-repmat(mean_all,1,size(data1,2))).^2,1))+mean(sum((data2-repmat(mean_all,1,size(data2,2))).^2,1)));

r = between/(within1+within2);