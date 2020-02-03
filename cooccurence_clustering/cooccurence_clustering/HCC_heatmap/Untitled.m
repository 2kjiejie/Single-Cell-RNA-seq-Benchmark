data = randn(200,50);
data(1:10,1:25) = data(1:10,1:25) + 10;
data(11:20,26:end) = data(11:20,26:end) + 4;
data(21:30,21:40) = data(21:30,21:40) - 4;
cancer_type = [ones(1,20), ones(1,5)*2, ones(1,15)*3, ones(1,10)*4];

similarity_matrix = corrcoef(data);

subplot(4,1,1:3); imagesc(similarity_matrix)
subplot(4,1,4); imagesc(cancer_type)

random_reorder = randsample(1:50,50);
data2 = data(:,random_reorder);
cancer_type2 = cancer_type(random_reorder);


%%%%%%%%%%%%%
data2
cancer_type2
similarity_matrix2 = corrcoef(data2);
[perm_r,perm_c] = HCC_heatmap(similarity_matrix2)
figure(2)
subplot(4,1,1:3); imagesc(similarity_matrix2(perm_r,perm_r))
subplot(4,1,4); imagesc(cancer_type2(perm_r))

