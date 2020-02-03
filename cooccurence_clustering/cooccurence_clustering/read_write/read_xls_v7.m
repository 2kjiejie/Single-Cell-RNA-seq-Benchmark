function [header,left_columns,data] = read_xls_v7(filename, sep_col, sep_row, skip_rows, left_column_ind, data_column_ind, num_header_rows, num_data_rows)  
% [header,left_columns,data] = read_xls_v7(filename,sep_col,sep_row,skip_rows,left_column_ind,data_column_ind,num_header_rows, num_data_rows)  
% restrictions: 
%   perfect rectangle matrix
% eg: 
%   [header,left_columns,data] = read_xls_v7(filename,char(9),char(10),0, [1,2], [], 2,inf)  


if ~exist('num_data_rows') || isempty('num_data_rows')
    num_data_rows = inf;
end
header = cell(0); left_columns = cell(0); data = [];
chunk_size = 50000;
fid = fopen(filename);

% figure out how many rows are there in total
total_num_rows = 0; d = ['aaaaa'];
while ~isempty(d)
    d = fread(fid,[1,chunk_size],'char'); d = char(d);
    num_complete_rows_in_chunk = sum(d==sep_row);
    total_num_rows = total_num_rows + num_complete_rows_in_chunk;
    if length(d)<chunk_size && length(d)~=0 && ~isequal(d(end),sep_row)
        total_num_rows = total_num_rows + 1;
    end
end
if total_num_rows==0, return; end
fseek(fid,0,-1);

%seek to right after the skipped rows
while skip_rows>0
    d = fread(fid,[1,chunk_size],'char'); d = char(d);
    if isempty(d) % trying to skip too many rows, more than what the file has
        fclose(fid);return 
    end
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk<1 % one chunk is too small, shorter than 2 lines
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
if ~exist('data_column_ind') || isempty(data_column_ind)
    data_column_ind = setdiff(1:num_columns, left_column_ind);
end



% begin reading the header rows
while num_header_rows>0 && ~isempty(d)
    d = fread(fid,[1,chunk_size],'char'); d = char(d); % read one chunk
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk<1  % if this chunk is too small 
        fseek(fid,-length(d),0); 
        chunk_size = chunk_size * 2;
        continue;
    end
    if num_complete_rows_in_chunk<num_header_rows
        ind = find(d==sep_row,1,'last');
        fseek(fid,-(length(d)-ind),0); % seek to right after the last complete row
    else
        ind = find(d==sep_row);
        ind = ind(num_header_rows);
        fseek(fid,-(length(d)-ind),0); % seek to right after the last header row
        num_complete_rows_in_chunk = num_header_rows;
    end        
    d(d==sep_row) = sep_col;
    tmp = regexp(d(1:(ind-1)),sep_col,'split'); 
    tmp2 = reshape(tmp,num_columns,num_complete_rows_in_chunk)';
    tmp3 = tmp2(:,[left_column_ind(:)', data_column_ind(:)']);
    header = [header;tmp3];
    num_header_rows = num_header_rows - num_complete_rows_in_chunk;
end


% begin reading the left_columns and data
if num_data_rows>total_num_rows-skip_rows-num_header_rows
    num_data_rows = total_num_rows-skip_rows-size(header,1);
end
data = zeros(num_data_rows,length(data_column_ind));
counter = 0;
fprintf('number of data rows remaining: %8d',num_data_rows);
while num_data_rows>0 && ~isempty(d)
    d = fread(fid,[1,chunk_size],'char'); d = char(d); % read one chunk
    if length(d)<chunk_size && ~isequal(d(end),sep_row)
        d = [d,sep_row];
    end
    num_complete_rows_in_chunk = sum(d==sep_row);
    if num_complete_rows_in_chunk<1  % if this chunk is too small 
        fseek(fid,-length(d),0); 
        chunk_size = chunk_size * 2;
        continue;
    end
    if num_complete_rows_in_chunk<num_data_rows
        ind = find(d==sep_row,1,'last');
        fseek(fid,-(length(d)-ind),0); % seek to right after the last complete row
    else
        ind = find(d==sep_row);
        ind = ind(num_data_rows);
        fseek(fid,-(length(d)-ind),0); % seek to right after the last header row
        num_complete_rows_in_chunk = num_data_rows;
    end        
    d(d==sep_row) = sep_col; 
    tmp = regexp(d(1:(ind-1)),sep_col,'split'); 
    tmp2 = reshape(tmp,num_columns,num_complete_rows_in_chunk)';
    left_columns = [left_columns; tmp2(:,left_column_ind)];
    s = cell2mat(strcat(reshape(tmp2(:,data_column_ind),1,length(data_column_ind)*num_complete_rows_in_chunk),{char(9)}));
    s(s==char(10))=[]; s(s==char(13))=[];
    tmp3 = reshape(str2num(s),num_complete_rows_in_chunk,length(data_column_ind));
    data(counter+1:counter+num_complete_rows_in_chunk,:) = tmp3;
    counter = counter + num_complete_rows_in_chunk;
    num_data_rows = num_data_rows - num_complete_rows_in_chunk;
    fprintf('\b\b\b\b\b\b\b\b%8d',num_data_rows);
end
fprintf('\n');

fclose(fid);
    