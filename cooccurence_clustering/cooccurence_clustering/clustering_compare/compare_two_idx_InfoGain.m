function dist = compare_two_idx_InfoGain(idx1,idx2)

idx1 = idx1(:);
idx2 = idx2(:);
idx = [idx1,idx2];
[B,I,J] = unique(idx,'rows'); % J is another idx that describes the "intersection" of the two idx's


info_gain_ratio1 = get_info_gain_ratio(J,idx1,2);
info_gain_ratio2 = get_info_gain_ratio(J,idx2,2);

dist = (nanmean(info_gain_ratio1) + nanmean(info_gain_ratio2))/2;
return


%%%%%%%%%%%%%%%%%%
function info_gain_ratio = get_info_gain_ratio(J,idx, min_allowable_size)

possible_idx = unique(idx);
info_gain_ratio = zeros(1, length(possible_idx));
for i=1:length(possible_idx)
    ind = find(idx==possible_idx(i));
    if length(ind)<min_allowable_size
        info_gain_ratio(i) = NaN;
    else
        info_gain_ratio(i) = entropy_idx(J(ind))/entropy_idx(1:min(max(J),length(ind)));
    end
end
