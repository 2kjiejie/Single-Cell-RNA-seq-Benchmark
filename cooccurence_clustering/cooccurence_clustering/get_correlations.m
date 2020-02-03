function r = get_correlations(all_genedata, class_labels,test_mode)
% r = get_correlations(all_genedata, class_labels,test_mode)
% can only handle two classes


if length(unique(class_labels))>2
    disp('P: sorry, this code can only handle two class!');
    return
end

% normalize to 0-mean 1-std
% all_genedata = all_genedata - repmat(mean(all_genedata,2),1,size(all_genedata,2));
% all_genedata = all_genedata./repmat(std(all_genedata')',1,size(all_genedata,2));
class_labels = class_labels ==max(class_labels); % transform class labels to 0,1
r = zeros(size(all_genedata,1),1);

if ~exist('test_mode')
    test_mode = 'correlation';
end

if ~isempty(intersect({test_mode},{'correlation'}))
    class_labels = class_labels - mean(class_labels);
    r = all_genedata*(class_labels');
    for i=1:length(r)
        if std(all_genedata(i,:))==0 
            r(i)=0;
        else
            r(i) = r(i)/(norm(all_genedata(i,:))*norm(class_labels));
        end
    end
    return
end

if ~isempty(intersect({test_mode},{'snr'}))
    data1 = all_genedata(:,class_labels==0);
    data2 = all_genedata(:,class_labels==1);
    mean_diff =  (mean(data2,2) - mean(data1,2));
    std1 = std(data1')'; 
    std2 = std(data2')'; 
    std_sum = std1 + std2;
    std_all = std(all_genedata')';
    std_sum(std_sum==0) = std_all(std_sum==0);    
    r = mean_diff./std_sum;
    ind = find(std(all_genedata')'==0); % channels that do not vary at all
    r(ind)=0;
    r(isinf(r))=0; 
    r(isnan(r))=0;
    return
end

if ~isempty(intersect({test_mode},{'ttest'}))
    data1 = all_genedata(:,class_labels==0);
    data2 = all_genedata(:,class_labels==1);
    n1 = size(data1,2);
    n2 = size(data2,2);
    r = (mean(data2,2) - mean(data1,2))./(sqrt(((n1-1)*var(data1')'+(n2-1)*var(data2')')/(n1+n2-2))*sqrt(1/n1+1/n2));
    ind = find(std(all_genedata')'==0); % channels that do not vary at all
    r(ind)=0;
    return
end
