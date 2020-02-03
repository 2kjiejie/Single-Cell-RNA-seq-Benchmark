function perm_data = randomly_shuffle_matrix(data)
% randomly shuffle rows of the data, shuffle each row independently

perm_data = zeros(size(data));
[~,I] = sort(rand(size(data)),2);
perm_data = data((I-1)*size(data,1)+repmat((1:size(data,1))',1,size(data,2)));
