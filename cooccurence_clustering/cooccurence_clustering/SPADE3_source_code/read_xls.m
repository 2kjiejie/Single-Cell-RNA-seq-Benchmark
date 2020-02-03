function table_content = read_xls(filename, sep, skip_rows, column_ind, num_rows)
% table_content = read_xls(filename,sep, skip_rows, column_ind, num_rows)

% sep=char(9);

% find out file size
fid = fopen(filename);
for i=1:skip_rows, tmp = fgetl(fid);end   % skip rows
tmp = fgetl(fid); line_count = 1; %content(line_count) = {tmp};
if tmp(end)~=sep && tmp(end-1)~=sep
    tmp=[tmp,sep];
end
separation_points = [0,find(tmp==sep)];
while ischar(tmp) 
    if isempty(tmp) || length(tmp)<length(separation_points)-1
        tmp = fgetl(fid); 
        continue;
    end
    tmp = fgetl(fid); line_count = line_count+1 %content(line_count) = {tmp};
    if exist('num_rows') && line_count>num_rows
        break;
    end
end
fclose(fid);
line_count = line_count-1;

if ~exist('column_ind') || isempty(column_ind)
    table_content = cell(line_count,length(separation_points)-1);
    column_ind = 1:size(table_content,2);
else
    table_content = cell(line_count,length(column_ind));
end

if exist('num_rows')
    table_content = table_content(1:num_rows,:);
end

fid = fopen(filename);
for i=1:skip_rows, tmp = fgetl(fid);end   % skip rows
for i=1:size(table_content,1)
    i
    tmp = fgetl(fid);
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
%     count = 1;
    for j=1:length(separation_points)-1
        if sum(column_ind==j)==1
            seg = tmp(separation_points(j)+1:separation_points(j+1)-1);
            count = find(column_ind==j);
            if isanumber(seg)==1 
                table_content(i,count) = {str2num(seg)};
            else
                table_content(i,count) = {seg};
            end
%             count = count + 1;
        end
    end    
end
fclose(fid);
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
if ~isempty(intersect(str,'+-')) && ~isempty(setdiff(find(str=='+' | str=='-'),1)) && ~isempty(setdiff(str(setdiff(find(str=='+' | str=='-'),1)-1),'eE'))
    result=0; return
end
if length(str2num(str))~=1
    result = 0; return
end
result = 1;
return