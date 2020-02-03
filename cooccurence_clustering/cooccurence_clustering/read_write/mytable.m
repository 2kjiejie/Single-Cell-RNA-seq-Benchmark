function [B, counts, table_stat] = table(data)
% function [B, counts, table_stat] = table(data)

data = data(:);

[B,I,J] = unique(data);

table_stat = cell(2,length(B));

for i=1:size(table_stat,2)
    if iscell(B(i))
        table_stat(1,i) = B(i);
    else
        table_stat{1,i} = B(i);
    end
    table_stat{2,i} = sum(J==i);
    counts(i) = sum(J==i);
end
return
