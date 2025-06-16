function [r2, ia, nse, kge, PBIAS, rmse, mae] = fitness_metrics(observed,modeled)
data1 = observed;
data2 = modeled;

% Calculate IA
numerator = sum((data1(:) - data2(:)).^2);
denominator = sum((abs(data2(:) - mean(data1(:))) + abs(data1(:) - mean(data1(:)))).^2);
ia = 1 - numerator / denominator;

% Calculate NSE
mean1 = mean(data1(:));
mean2 = mean(data2(:));
nse = 1 - sum((data1(:) - data2(:)).^2) / sum((data1(:) - mean1).^2);

% PBIAS
PBIAS = sum(data1(:) - data2(:))/sum(data1(:)); 

% Calculate KGE
std1 = std(data1(:));
std2 = std(data2(:));
r = corrcoef(data1(:), data2(:));
r = r(1, 2);
alpha = std2 / std1;
beta = mean2/mean1;
kge = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);


% Calculate RMSE
rmse = sqrt(nanmean((data2 - data1).^2));

% Calculate MAE
mae = nanmean(abs(data2 - data1));

% Calculate R²
r2 = corrcoef(data1, data2, 'Rows','complete');
r2 = r2(1, 2)^2; % Extract R² value

metric_values = [r2, ia, nse, kge, PBIAS, rmse, mae];
end