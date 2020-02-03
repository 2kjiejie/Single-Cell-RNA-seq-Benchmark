function [header, feature_names, data] = read_xls_v4(filename) 
% [header, feature_names, data] = read_xls_v4(filename) 
% restrictions: 
%  (1) one row of header
%  (2) one column of left_column
%  (3) sep == char(9)

    fid = fopen(filename);
    d = fread(fid,[1,inf],'char'); d = char(d); fclose(fid);
    d(d==char(13))=[];
    % header
    header = regexp(d(1:find(d==char(10),1)-1),char(9),'split');
    num_columns = length(regexp(d(1:find(d==char(10),1)),char(9),'split'));
    % feature names
    ind1 = find(d==char(10)); ind1(end)=[]; ind1 = [0,ind1];
    ind2 = find(d==char(9)); ind2 = ind2(1:(num_columns-1):end);
    indicator = repmat(logical(0),1,length(d));
    for i=1:length(ind1)
        indicator(ind1(i)+1:ind2(i))=1;
    end
    feature_names = d(indicator==1); feature_names(end)=[];
    feature_names = regexp(feature_names,char(9),'split')'; %including the first element of the header row
    feature_names(1)=[]; %exclude the first element of the header row
    clear ind1 ind2 indicator
    % data
    ind1 = find(d==char(9)); ind1 = ind1(1:(num_columns-1):end); ind1(1)=[];
    ind2 = find(d==char(10));ind2(1)=[];
    data = zeros(length(ind1),num_columns-1);
    for i=1:length(ind1)
        [i,length(ind1)]
        data(i,:) = str2num(d(ind1(i)+1:ind2(i)-1));
    end
    clear d
    
    