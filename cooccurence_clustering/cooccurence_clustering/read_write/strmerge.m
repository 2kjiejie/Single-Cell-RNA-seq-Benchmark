function [str] = strmerge(seg, delimiter)
% [str] = strmerge(seg, delimiter)
%       seg is a vector of cell strings
%       delimiters can be one string, or a char

str=[];
for i=1:length(seg)
    str = [str, seg{i}];
    if i~=length(seg)
        str = [str, delimiter];
    end
end
