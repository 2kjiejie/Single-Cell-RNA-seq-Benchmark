function [score] = compare_ranking_lists(I1,I2,sort_mode)
% sort_mode = 'ascend' or 'descend'


if sum(abs(I1-round(I1)))~=0 &&  sum(abs(I2-round(I2)))~=0 
    if ~exist('sort_mode') || isempty(intersect({'ascend','descend'},sort_mode))
        sort_mode = 'descend';
    end
    [Y,I1] = sort(I1,sort_mode);
    [Y,I2] = sort(I2,sort_mode);
end
    

I1 = I1(:); I2 = I2(:);
[Y,I3] = sort(rand(size(I1)),'descend'); 
[Y,I4] = sort(rand(size(I1)),'descend'); 
score=0;null_score = 0;
x = 1:ceil(length(I1)/100):length(I1);
for i=1:length(x)
    o(i) = length(intersect(I1(1:x(i)),I2(1:x(i))));
    n(i) = length(intersect(I3(1:x(i)),I4(1:x(i))));
    score = score + (1-o(i)/x(i));
    null_score = null_score + (1-n(i)/x(i));
end
score = score/length(x);
null_score = null_score/length(x);
score = [score,null_score];
plot([1,length(I1)],[1,length(I1)],'g');hold on;
plot(x,o,'b',x,n,'g')
hold off