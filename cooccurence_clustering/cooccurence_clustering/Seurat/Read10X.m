function [data, gene_names, cell_names] = Read10X(folder_name)

% filenames = getfilenames(folder_name);

% cell_names = read_xls_v2(fullfile(folder_name,'barcodes.tsv'),char(9),0);
% gene_names = read_xls_v2(fullfile(folder_name,'genes.tsv'),char(9),0);
% gene_names = gene_names(:,2);
% sparse_data = read_xls_v2(fullfile(folder_name,'matrix.mtx'),' ',2);

% cell names
fprintf('Loading cell names ... ');
tic
fid = fopen(fullfile(folder_name,'barcodes.tsv'));
d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
d(d==char(13))=[];
cell_names = regexp(d,char(10),'split');
cell_names(end)=[];
toc

% gene names
fprintf('Loading gene names ... ');
tic
fid = fopen(fullfile(folder_name,'genes.tsv'));
d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
d(d==char(13))=[];
d(d==char(9)) = char(10);
gene_names = regexp(d,char(10),'split');
gene_names = gene_names(2:2:end)';
toc

% raw_data
fprintf('Loading raw data ... ');
tic
fid = fopen(fullfile(folder_name,'matrix.mtx'));
d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
d(1:max(find(d==char(10),2)))=[];
d(d==char(10))=char(9);
s = str2num(d);
num_genes = s(1);
num_cells = s(2);
num_entries = s(3);
raw_data = reshape(s(4:end),3,length(s)/3-1);
toc

% reformat
fprintf('Converting data of %d genes * %d cells into matrix ... ',num_genes, num_cells);
tic
data = sparse(raw_data(1,:),raw_data(2,:), raw_data(3,:));
if size(data,1)~=num_genes || size(data,2)~=num_cells
    data(num_genes,num_cells)=0;
end
toc