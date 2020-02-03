function [outfile] = write_symmetric_adj_into_graph_file(A, outfile)

if ~exist('outfile')
    outfile = 'G.txt';
end

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