function [data, sample_names, gene_names] = read_expr(filename, sep, sample_row, data_start_row, gene_column, data_column, num_rows)

% this program is only suitable for reading one column of gene discription

% sep=char(9);

% find out file size
fid = fopen(filename);

if isempty(sample_row) || sample_row ==0
    sample_names = [];
else
    skip_rows = sample_row - 1;
    if skip_rows~=0
        for i=1:skip_rows, tmp = fgetl(fid);end   % skip rows
    end
    tmp = fgetl(fid); % this row contains sample names
    if tmp(end)~=sep && tmp(end-1)~=sep
        tmp=[tmp,sep];
    end
    separation_points = [0,find(tmp==sep)];
    for j=1:length(separation_points)-1
        seg = tmp(separation_points(j)+1:separation_points(j+1)-1);
        sample_names{1,j} = seg;
    end
    sample_names = sample_names(data_column);
end

% start to read the data
if exist('num_rows') && ~isempty(num_rows)
    data = zeros(num_rows,length(data_column));
else
    data = zeros(10000,length(data_column));
end
gene_names = [];
skip_rows = data_start_row - sample_row - 1; % there might be rows between sample names and numerical data
if skip_rows>0
    for i=1:skip_rows, tmp = fgetl(fid);end   % skip rows
end
line_count = 0;
fprintf('reading feature ... '); fprintf('%8d%',0);
tmp = 'ch'; % ch doesn't mean anything, just make sure tmp is char
while ischar(tmp) 
    tmp = fgetl(fid);
    if ~ischar(tmp) && isnumeric(tmp) && tmp==-1
        break;
    end
    if length(tmp)>=24 && sum(tmp(1:24)=='!series_matrix_table_end')==24
        table_content = table_content(1:i-1,:);
        break;
    end
    while isempty(tmp) || length(tmp)<length(separation_points)-1
        tmp = fgetl(fid);
    end
    if tmp(end)~=sep
        tmp=[tmp,sep];
    end
    separation_points = [0,find(tmp==sep)];
    seg = tmp(separation_points(gene_column)+1:separation_points(gene_column+1)-1);
    gene_names = [gene_names; {seg}];
    data_tmp = zeros(1,max(data_column));
    for j=(data_column(:)')
        if separation_points(j) - separation_points(j+1) == -1
            data_tmp(j) = NaN;
        else
            data_tmp(j) = str2num(tmp(separation_points(j)+1:separation_points(j+1)-1));
        end
    end
    line_count = line_count + 1;
    data(line_count,:) = data_tmp(data_column);
    fprintf('\b\b\b\b\b\b\b\b%8d%',line_count);
    if exist('num_rows') && ~isempty(num_rows) && num_rows<=line_count
        break;
    end
end
fprintf('\n');
fclose(fid);
data = data(1:line_count,:);
return


function [result] = isanumber(str)
if isempty(str)
    result=0; return
end
while (str(1)==char(9) || str(1)==32 || str(1)==13 || str(1)==',') && length(str)>=2 % 9--tab, 32-space, 13-enter
    str = str(2:end);
end
while (str(end)==char(9) || str(end)==32 || str(end)==13 || str(end)==',') && length(str)>=2 % 9--tab, 32-space, 13-enter
    str = str(1:end-1);
end
if ~isempty(setdiff(str, '0123456789.+-E'))
    result = 0; return
end
if length(str2num(str))~=1
    result = 0; return
end
result = 1;
return