function [seg] = strsplit(str, delimiters)
% [seg] = strsplit(str, delimiters)
%       delimiters can be one string, or a cell vector of strings
%       if one delimiter is a substr of another, the result might come out weird


all_ind = [];
delimiter_len = [];
if iscell(delimiters)
    for i=1:length(delimiters)
        ind = strfind(str,delimiters{i});
        all_ind = [all_ind,ind];
        delimiter_len = [delimiter_len, length(delimiters{i})*ones(1,length(ind))];
    end
else
    ind = strfind(str,delimiters);
    all_ind = ind;
    delimiter_len = [delimiter_len, length(delimiters)*ones(1,length(ind))];
end
[Y,I] = sort(all_ind);
all_ind = all_ind(I);
delimiter_len = delimiter_len(I);

segment_start = [1,all_ind+delimiter_len];
segment_end = [all_ind-1,length(str)];

seg=[];
for i=1:length(segment_start)
    seg{i} = str(segment_start(i):segment_end(i));
end

return

