function [table_content] = read_xls_v6(filename,sep_col,sep_row,skip_rows,column_ind,num_rows_to_read) 
% [table_content] = read_xls_v6(filename,sep_col,sep_row,skip_rows,column_ind,num_rows_to_read)  
% restrictions: 
%   perfect rectangle matrix
% eg: 
%   [table_content] = read_xls_v6(filename,char(9),char(10),0,[],inf)  


if ~exist('num_rows_to_read') || isempty('num_rows_to_read')
    num_rows_to_read = inf;
end
table_content=cell(0);
chunk_size = 50000000;
fid = fopen(filename);

%seek to right after the skipped rows
while skip_rows>0
    d = fread(fid,[1,chunk_size],'char'); d = char(d);
    if isempty(d) % trying to skip too many rows, more than what the file has
        fclose(fid);return 
    end
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk<=1 % one chunk is too small, shorter than 2 lines
        fseek(fid,-length(d),0); % reverse the fread a few lines above
        chunk_size = chunk_size * 2;
        continue;
    end
    if num_complete_rows_in_chunk<skip_rows
        fseek(fid,-(length(d)-find(d==sep_row,1,'last')),0); % seek to right after the last complete row
        skip_rows = skip_rows - num_complete_rows_in_chunk;
    else % meaning that num_complete_rows_in_chunk>=skip_rows
        ind = find(d==sep_row);
        fseek(fid,-(length(d)-ind(skip_rows)),0); % seek to right after the last skipped row
        skip_rows = 0;
    end
end


% read the first line
num_complete_rows_in_chunk=0;
while num_complete_rows_in_chunk==0
    d = fread(fid,[1,chunk_size],'char'); d = char(d);
    fseek(fid,-length(d),0); 
    if isempty(d)
        fclose(fid); return
    end
    if length(d)<chunk_size && ~isequal(d(end),sep_row)
        d = [d,sep_row];
    end
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk==0 
        chunk_size = chunk_size * 2;
    end
end
num_columns = sum(d(1:(find(d==sep_row,1)-1))==sep_col)+1;
if ~exist('column_ind') || isempty(column_ind)
    column_ind = 1:num_columns;
end



% begin reading the data
while num_rows_to_read>0 && ~isempty(d)
    d = fread(fid,[1,chunk_size],'char'); d = char(d); % read one chunk
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk<=1  % if this chunk is too small 
        fseek(fid,-length(d),0); 
        chunk_size = chunk_size * 2;
        continue;
    end
    if num_complete_rows_in_chunk<num_rows_to_read
        ind = find(d==sep_row,1,'last');
        fseek(fid,-(length(d)-ind),0); % seek to right after the last complete row
    else
        ind = find(d==sep_row);
        ind = ind(num_rows_to_read);
        fseek(fid,-(length(d)-ind),0); % seek to right after the last complete row
        num_complete_rows_in_chunk = num_rows_to_read;
    end        
%     d = d(1:(ind-1));
%     d(d==sep_row) = sep_col;
%     tmp = regexp(d,sep_col,'split'); 
%     tmp = reshape(tmp,num_columns,num_complete_rows_in_chunk)';
%     tmp = tmp(:,column_ind);

    d(d==sep_row) = sep_col;
    tmp = regexp(d(1:(ind-1)),sep_col,'split'); 
    tmp2 = reshape(tmp,num_columns,num_complete_rows_in_chunk)';
    tmp3 = tmp2(:,column_ind);
    table_content = [table_content;tmp3];
    num_rows_to_read = num_rows_to_read - num_complete_rows_in_chunk;
end

fclose(fid);
    