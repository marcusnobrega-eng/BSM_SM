function [lags, cross_corr] = correlation_analysis(x, y, max_lag)
    % Function to compute correlation and lag correlation between two variables
    % 
    % Inputs:
    %   x       - First time series (vector)
    %   y       - Second time series (vector)
    %   max_lag - Maximum lag for cross-correlation
    %
    % Outputs:
    %   - Plots correlation results

    % Ensure both variables are column vectors
    x = x(:);
    y = y(:);

    % Compute cross-correlation (lag correlation)
    [cross_corr, lags] = xcorr(x - mean(x), y - mean(y), max_lag, 'coeff');

end
