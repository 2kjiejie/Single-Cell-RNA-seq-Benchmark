function table_content = StrCellCat(table_content,dim,insert_str)

possible_new_delimiter = {'`','!','@','$','&','*','(',')','/','\'};

% control for unreasonable input 
if length(size(table_content))>2
    error('dimensionality of table_content is too big, more than 2');    
end

% figure out which delimiator is useable, meaning not exist in the table_content variable  
tmp = cell2mat(table_content(:)');
new_delimiter=[];
for i=1:length(possible_new_delimiter)
    if isempty(strfind(tmp, possible_new_delimiter{i})) && isempty(strfind(insert_str, possible_new_delimiter{i}))
        new_delimiter = possible_new_delimiter{i};
        break;
    end
end
if isempty(new_delimiter)
    error('table_content is too complex, containing all possible new delimiters');
end


% two case of cat vertically or horizontally
if dim == 1
    tmp = table_content;
    if size(table_content,1)==1
        return
    end
%     tmp(1:end-1,:) = strcat(tmp(1:end-1,:),insert_str);
%     tmp(end,:) = strcat(tmp(end,:),new_delimiter);
%     tmp = cell2mat(tmp(:)');
%     tmp = regexp(tmp,new_delimiter,'split');
%     tmp(end) = [];
%     table_content = tmp;
    tmp2 = cell(size(tmp,1)*2, size(tmp,2));
    tmp2(1:2:end,:) = tmp;
    tmp2(2:2:end-2,:) = repmat({insert_str},size(tmp,1)-1, size(tmp,2));
    tmp2(end,:) = repmat({new_delimiter}, 1, size(tmp,2));
    tmp2 = cell2mat(tmp2(:)');
    tmp2 = regexp(tmp2,new_delimiter,'split');
    tmp2(end) = [];
    table_content = tmp2;   
elseif dim==2
    tmp = table_content';
    if size(table_content,2)==1
        return
    end
%     tmp(1:end-1,:) = strcat(tmp(1:end-1,:),insert_str);
%     tmp(end,:) = strcat(tmp(end,:),new_delimiter);
%     tmp = cell2mat(tmp(:)');
%     tmp = regexp(tmp,new_delimiter,'split');
%     tmp(end) = [];
%     table_content = tmp';    
    tmp2 = cell(size(tmp,1)*2, size(tmp,2));
    tmp2(1:2:end,:) = tmp;
    tmp2(2:2:end-2,:) = repmat({insert_str},size(tmp,1)-1, size(tmp,2));
    tmp2(end,:) = repmat({new_delimiter}, 1, size(tmp,2));
    tmp2 = cell2mat(tmp2(:)');
    tmp2 = regexp(tmp2,new_delimiter,'split');
    tmp2(end) = [];
    table_content = tmp2';   
else
    error('dim can only be 1 or 2');
end
    
    


