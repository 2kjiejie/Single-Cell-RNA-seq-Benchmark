function r = compute_snr(data1, data2)

% each row is detection rate of one pathway
% each column is one cell

mean1 = mean(data1,2);
mean2 = mean(data2,2);
mean_diff =  (mean(data2,2) - mean(data1,2));

std1 = std(data1')'; 
std2 = std(data2')'; 
std_all = std([data1,data2]')';

std_sum = std1 + std2;
std_sum(std_sum==0) = std_all(std_sum==0);    

std_min = min(std1,std2);
std_min(std_min==0) = std_sum(std_min==0);

r = abs(mean_diff)./(std_min*2);

ind = find(std_all==0); % channels that do not vary at all
r(ind)=0;
r(isinf(r))=0; 
r(isnan(r))=0;
