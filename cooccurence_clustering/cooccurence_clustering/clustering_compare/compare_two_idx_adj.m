function [a,b] = compare_two_idx_adj(idx1,idx2)

idx1 = idx1(:);
adj1 = repmat(logical(0), length(idx1), length(idx1));
for i=1:length(idx1)
    adj1(:,i) = idx1 == idx1(i);
end
% adj1 = (repmat(idx1,1,length(idx1)) == repmat(idx1',length(idx1),1));


idx2 = idx2(:);
adj2 = repmat(logical(0), length(idx2), length(idx2));
for i=1:length(idx2)
    adj2(:,i) = idx2 == idx2(i);
end
% adj2 = (repmat(idx2,1,length(idx2)) == repmat(idx2',length(idx2),1));


a = sum(sum(adj1==adj2))/prod(size(adj1));
b = sum(sum(adj1 & adj2 ==1))/sum(sum(adj1==1 | adj2==1));

