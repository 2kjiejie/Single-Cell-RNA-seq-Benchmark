function [A, graph_file]=Seurat_SNN(d, k, snn_prune, distance_metric, outfile)
% data input. row: object/sample; column: attribute.
% outfile: strength output
% k: nearest neighbor
% distance: methods for calculation distance (euclidean, correlation ...)

if ~exist('k')
    k=30;
end

if ~exist('snn_prune')
    snn_prune = 1/15;
end

if ~exist('distance_metric')
    distance_metric='euclidean';
end

if ~exist('outfile')
    outfile = 'G.txt';
end

numSpl=size(d,1);
tic
IDX=knnsearch(d, d, 'K', k, 'Distance', distance_metric);
% fprintf('Constructing SNN graph ... %3d%%', 0);
% directed_adj=zeros(numSpl,numSpl);
% for i=1:numSpl
%     directed_adj(i,IDX(i,:))=1; 
%     fprintf('\b\b\b\b%3d%%', round(i/numSpl*100));
% end
fprintf('Constructing SNN graph ... ');
directed_adj = sparse(reshape(repmat(IDX(:,1),1,size(IDX,2)),prod(size(IDX)),1), IDX(:), ones(prod(size(IDX)),1));
A = directed_adj*directed_adj';
tmp_ind = find(A~=0);
A(tmp_ind) = A(tmp_ind)./(k+(k-A(tmp_ind)));
tmp_ind = tmp_ind(A(tmp_ind)<snn_prune);
A(tmp_ind)=0;
toc

tic
fprintf('Writing graph into file ... %3d%%', 0);
fid=fopen(outfile,'w');
B = tril(A,-1);
[a,b] = find(B~=0);
tmp = [b,a,B(B~=0)];
[~,I] = sort(b);
tmp = tmp(I,:); 
tmp(:,1:2) = tmp(:,1:2)-1;
block_size = 1000;
sep = char(9);
for i=1:block_size:size(tmp,1)
    tmp_block = tmp(i:min(i+1000,end),:);
    data_tmp_str = num2str(reshape(tmp_block',1,prod(size(tmp_block))));
    data_tmp_str(regexp(num2str(data_tmp_str),' +','start'))=sep;
    if data_tmp_str(1)==sep
        data_tmp_str(1)=[];
    end
    data_tmp_str(data_tmp_str==' ') = [];
    ind = find(data_tmp_str==sep);
    data_tmp_str(ind(3:3:end))=char(10);
    fwrite(fid,data_tmp_str,'char');
    fprintf('\b\b\b\b%3d%%', round((i+size(tmp_block,1))/size(tmp,1)*100));
end
fclose(fid);
toc

graph_file = outfile;

return
