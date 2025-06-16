function [metrics] = calculate_metrics(observed, simulated)
    % calculate_metrics: Calculates NSE, RMSE, PBIAS, R², and KGE for observed and simulated data.
    %
    % Inputs:
    %   observed  - vector of observed values
    %   simulated - vector of simulated or predicted values
    %
    % Output:
    %   metrics - a struct containing NSE, RMSE, PBIAS, R², and KGE
    
    % Ensure inputs are column vectors
    observed = observed(:);
    simulated = simulated(:);
    
    % Error checks
    if length(observed) ~= length(simulated)
        error('Observed and simulated vectors must have the same length.');
    end
    
    if any(isnan(observed)) || any(isnan(simulated))
        error('Input data contains NaN values. Please clean the data.');
    end
    
    % Mean values
    obs_mean = mean(observed);
    sim_mean = mean(simulated);
    
    % 1. Nash-Sutcliffe Efficiency (NSE)
    numerator = sum((observed - simulated).^2);
    denominator = sum((observed - obs_mean).^2);
    NSE = 1 - (numerator / denominator);
    
    % 2. Root Mean Square Error (RMSE)
    RMSE = sqrt(mean((observed - simulated).^2));
    
    % 3. Percent Bias (PBIAS)
    PBIAS = 100 * sum(observed - simulated) / sum(observed);
    
    % 4. R-squared (R²)
    r = corr(observed, simulated); % Correlation coefficient
    R_squared = r^2;
    
    % Alternatively, calculate R² using correlation:
    % R_squared = (corr(observed, simulated))^2;
    
    % 5. Kling-Gupta Efficiency (KGE)
    % Components of KGE
    r = corr(observed, simulated); % Correlation coefficient
    alpha = std(simulated) / std(observed); % Variability ratio
    beta = sim_mean / obs_mean; % Bias ratio
    KGE = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
    
    % Output metrics as a struct
    metrics = struct(...
        'NSE', NSE, ...
        'RMSE', RMSE, ...
        'PBIAS', PBIAS, ...
        'R_squared', R_squared, ...
        'KGE', KGE);
end
