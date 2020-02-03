function ind = isInListFront(all_filenames,frontstring)
% function ind = isInList(all_filenames,frontstring)
% 

all_filenames = upper(all_filenames);
frontstring = upper(frontstring);
ind = [];
for i=1:length(all_filenames)
    if length(all_filenames{i})>=length(frontstring) && isequal(all_filenames{i}(1:length(frontstring)),frontstring)
       ind = [ind,i];
    end
end

