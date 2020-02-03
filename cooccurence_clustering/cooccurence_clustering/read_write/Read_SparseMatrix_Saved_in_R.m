function [data, row_names, col_names] = Read_SparseMatrix_Saved_in_R(filename_matrix, filename_row, filename_col)
%
% Assume the data files are generated from R (eg. the following R code)
%     library(tidyverse)
%     tm.droplet.matrix = readRDS("TM_droplet_mat.rds")
%     writeMM(tm.droplet.matrix,'TM_droplet_mat.txt')
%     write.csv(rownames(tm.droplet.matrix),'TM_droplet_rownames.csv')
%     write.csv(colnames(tm.droplet.matrix),'TM_droplet_colnames.csv')
% 
% Example use of this code: 
%     Read_SparseMatrix_Saved_in_R('TM_droplet_mat.txt', 'TM_droplet_rownames.csv', 'TM_droplet_colnames.csv')
% 

% row names
fprintf('Loading row names ... ');
tic
fid = fopen(filename_row);
d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
d(d==char(13))=[];
tmp = regexp(d,'"','split');
row_names = tmp(8:4:end);
toc

% col names
fprintf('Loading col names ... ');
tic
fid = fopen(filename_col);
d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
d(d==char(13))=[];
tmp = regexp(d,'"','split');
col_names = tmp(8:4:end);
toc


% raw_data
fprintf('Loading raw data ... ');
tic
% fid = fopen(filename_matrix);
% d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
% d(d==char(13))=[]; 
% d(1:max(find(d==char(10),1)))=[];
% 
% lines_end = find(d==char(10))-1;
% lines_start = [1,lines_end(1:end-1)+2];
% 
% s = d(lines_start(1):lines_end(1));
% num_rows = s(1);
% num_cols = s(2);
% num_entries = s(3);
% 
% block_size = 100000;
% raw_data = zeros(length(lines_start)-1,3);
% fprintf('%3d%%',0);
% for i=2:block_size:length(lines_start)
%     lines_to_convert = i:min(i+block_size-1,length(lines_start));
%     raw_data(lines_to_convert-1,:) = str2num(d(lines_start(lines_to_convert(1)):lines_end(lines_to_convert(end))));
%     fprintf('\b\b\b\b%3d%%',round(lines_to_convert(end)/length(lines_start)*100));
% end

s = importdata(filename_matrix);
num_rows = s.data(1,1);
num_cols = s.data(1,2);
num_entries = s.data(1,3);
raw_data = s.data(2:end,:);
clear s
toc

% reformat
fprintf('Converting data of %d rows * %d columns into matrix ... ',num_rows, num_cols);
tic
data = sparse(raw_data(:,1),raw_data(:,2), raw_data(:,3));
if size(data,1)~=num_rows || size(data,2)~=num_cols
    data(num_rows,num_cols)=0;
end
toc