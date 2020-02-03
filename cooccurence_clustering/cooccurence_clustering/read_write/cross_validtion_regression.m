function [y_hat, mae, y] = cross_validtion_regression(image_level_cluster_count, all_para, varargin)

rng(1)
% Handle Parameters
% K: number of clusters in k-means [default:10]
% Folds: number of folds in cross validation [default: 5]
% ynames: whether or not to show figure about the regression [default: true]
% log: whether or not to log trnasform the response [default: false]
pnames = { 'K', 'Folds', 'Showfigure', 'ynames', 'log'};
dflts  = { 10,5, false, {'gx', 'kb'}, false};
[K, Folds, showfigure, ynames, ylog] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if ischar(K)
    [K, status] = str2num(K);
    if ~status || round(K) ~= K
        msgID = 'MATLAB:ValueError';
        msgtext = 'The K you provided is not a valid integer';
        ME = MEception(msgID, msgtext);
        throw(ME)
    end
end
if ischar(Folds)
    [Folds, status] = str2num(Folds);
    if ~status || round(Folds) ~= Folds
        msgID = 'MATLAB:ValueError';
        msgtext = 'The Folds you provided is not a valid integer';
        ME = MEception(msgID, msgtext);
        throw(ME)
    end
end    
%% Prepare the data

x = image_level_cluster_count;

% Decide whether or not log transform the data
if ylog
    y1 = log(all_para(:,2)); % gx
    %y1 = (y1-min(y1))./range(y1);
    y2 = log(all_para(:,3)); % kb
    %y2 = (y2-min(y2))./range(y2);
else
    y1 = all_para(:,2); % gx
    y2 = all_para(:,3); % kb
    %y1 = all_para(:,2)./max(all_para(:,2)); % gx
    %y2 = all_para(:,3)./max(all_para(:,3)); % kb
end
%% Cross Validation
y = [y1 y2];
y_hat = zeros(size(y));
indices = crossvalind('Kfold', size(y_hat, 1) ,Folds);
for i = 1:Folds
    cvidx = (indices==i);
    train_x = [ones(sum(~cvidx),1) x(~cvidx,:)];
    train_y = y(~cvidx, :);
    test_x = [ones(sum(cvidx),1) x(cvidx,:)];
    [beta, ~] = mvregress(train_x, train_y);
    y_hat(cvidx,:) = test_x * beta;
end
if showfigure
    % show 1D difference
    figure 
    subplot(1,2,1)
    plot(1:length(y(:,2)), y(:,2), 1:length(y(:,2)), y_hat(:,2));
    ylabel(sprintf('normalized %s', ynames{2}))
    legend({'Truth','Predicted'},'Location','best')
    xlabel('sample')
    subplot(1,2,2)
    plot(1:length(y(:,1)), y(:,1), 1:length(y(:,1)), y_hat(:,1));
    ylabel(sprintf('normalized %s', ynames{1}))
    legend({'Truth','Predicted'},'Location','best')
    xlabel('sample')
    
    % show how different in 2D 
    figure 
    plot(y(:,1), y(:,2),'o', y_hat(:,1), y_hat(:,2),'o');
    legend({'Truth', 'Predicted'},'Location','best')
    hold on
    for i = 1:size(y,1)
        plot([y(i,1) y_hat(i,1)], [y(i,2) y_hat(i,2)], 'k-');
    end
    xlabel('gx')
    ylabel('kb')
end
% mae = median(abs(y_hat-y)./abs(y));
mae = mean(abs(y_hat-y));





