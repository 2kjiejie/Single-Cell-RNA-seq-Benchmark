function [table_content] = read_xls_v5(filename) 
% [table_content] = read_xls_v4(filename) 
% restrictions: 
%   sep == char(9)
%   perfect rectangle matrix

    fid = fopen(filename);
    d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
    d(d==char(13))=[];
    % first row and number of columns
    header = regexp(d(1:find(d==char(10),1)-1),char(9),'split');
    num_columns = length(regexp(d(1:find(d==char(10),1)),char(9),'split'));
    % number of rows
    num_rows = sum(d==char(10));
    % get the data
    d(d==char(10))=char(9);
    table_content = regexp(d,char(9),'split');
    table_content = reshape(table_content(1:end-1),num_columns,num_rows)';
    clear d
    
    