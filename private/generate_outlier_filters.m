function [ FILTERS ] = generate_outlier_filters(X)
% This function will z-score the row and column means and remove those
% which have a z-score greater than 5.  The process it performed on
% columns, and then rows.  Each time outliers are removed, the reduced data
% is checked again for outliers, until no outliers remain.

    FILTERS = struct(...
        'label',{'rowfilter','colfilter'},...
        'dimension',{1,2},...
        'filter', {[],[]}...
    );

    % COLUMNS
    z = abs(zscore(mean(X,1))) > 5;
    outliers = z;
    while nnz(z) > 0
        z = abs(zscore(mean(X(:,~outliers),1))) > 5;
        outliers(~outliers) = outliers(~outliers) | z;
    end
    FILTERS(2).filter = ~outliers;
    
    % ROWS
    z = abs(zscore(mean(X,2))) > 5;
    outliers = z;
    while nnz(z) > 0
        z = abs(zscore(mean(X(~outliers,:),2))) > 5;
        outliers(~outliers) = outliers(~outliers) | z;
    end
    FILTERS(1).filter = ~outliers;
end