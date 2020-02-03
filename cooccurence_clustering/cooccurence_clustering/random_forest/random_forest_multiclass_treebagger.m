function [accuracy, auc, confusion_table, cv_pred_Y, cv_pred_prob] = random_forest_multiclass_treebagger(X,Y,feature_names, isdisplay)

Y=Y(:)';
if size(X,2)==length(Y)
    X = X';
end

if ~exist('feature_names')
    feature_names = cell(1,size(X,2));
    for i=1:length(feature_names)
        feature_names{i} = ['x_',num2str(i)];
    end
end

if ~exist('isdisplay') || ~isequal(isdisplay,1)
    isdisplay=0;
end



ntree = 50; 
B = TreeBagger(ntree,X,Y,'InBagFraction',1,'Method','classification','OOBPrediction','off','OOBPredictorImportance','off');
[Yfit,scores,stdevs] = predict(B,X);
score_reorder = str2double(unique(Yfit));
preds(:,score_reorder) = scores;


cv_ind = crossvalind('Kfold', Y, 5);
preds2 = zeros(size(preds));
for i=1:5
    train_X = X(cv_ind~=i,:);
    train_Y = Y(cv_ind~=i);
    test_X = X(cv_ind==i,:);
    test_Y = Y(cv_ind==i);
    B2 = TreeBagger(ntree,train_X,train_Y(:),'InBagFraction',1,'Method','classification','OOBPrediction','off','OOBPredictorImportance','off');
    [Yfit,scores,stdevs] = predict(B2,test_X);
    preds2(cv_ind==i,score_reorder) = scores;
end
cv_pred_prob = preds2;



auc1 = [];
for i=1:max(Y)
    [~,~,~,tmp] = perfcurve(Y==i',preds(:,i)',1);
    auc1 = [auc1,tmp];
end

auc2 = [];
for i=1:max(Y)
    [~,~,~,tmp] = perfcurve(Y==i',preds2(:,i)',1);
    auc2 = [auc2,tmp];
end

[~,tmp]=max(preds,[],2);
accuracy1 = mean(tmp == Y(:));
[~,tmp]=max(preds2,[],2);
accuracy2 = mean(tmp == Y(:));
cv_pred_Y = tmp;

if isdisplay==1
    accuracy1
    accuracy2
    auc1
    auc2
end

accuracy = accuracy2;
auc = auc2;
confusion_table = confusionmat( double(tmp),double(Y));

