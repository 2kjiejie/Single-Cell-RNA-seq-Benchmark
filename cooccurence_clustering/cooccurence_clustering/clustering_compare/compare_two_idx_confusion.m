function [C] = compare_two_idx_confusion(idx1,idx2)

idx1 = standardize_idx(idx1(:)); % delete those unassigend (0s) in the indicator vector
idx2 = standardize_idx(idx2(:));
unique_idx1 = unique(idx1');
unique_idx2 = unique(idx2');

C = zeros(length(unique_idx1), length(unique_idx2));
for i=unique_idx1
    for j=unique_idx2
        C(i,j) = sum(idx1==unique_idx1(i) & idx2==unique_idx2(j)); 
    end
end
