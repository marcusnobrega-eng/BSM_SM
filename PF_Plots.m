%% Particle Filter Plots



flag_PF = 1;
if flag_PF == 1
% Particle Filter Data
dischargeStore_PF(dischargeStore_PF > 100) = nan; % These are values with wrong simulations stored
time_step = tSave(1,1);                 % Seconds
time_plot_days = tSave/(3600*24);       % days

close all
% load z_estimated.csv
obs_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Observed_Hydrograph.csv';
z_estimated = table2array(readtable(obs_dir));
discharge_data = z_estimated(:,2);
discharge_data = discharge_data(1:size(dischargeStore_PF,1));

bechmark_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Calibration\calibrated_discharge.csv';
z_estimated = table2array(readtable(bechmark_dir));
discharge_benchmark = z_estimated(:,2); % benchmark discharge

z_estimated = nanmean(dischargeStore_PF');
std_estimated_z = nanstd(dischargeStore_PF');

[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps();
figure(1);
subplot(8,1,[1,2])
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, .4, 1]);
plot(time_plot_days, discharge_data, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_z, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(time_plot_days, z_estimated, '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([time_plot_days, fliplr(time_plot_days)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
% xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); 
ylabel('Discharge [mm/h]','Interpreter','latex','FontSize',16)
plot(time_plot_days, discharge_benchmark, '-','linewidth',1.5,'Color',pallete.red_colors(2,:));
% xlim([7,8])
ylim([0,2]);
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Observed','Particle filter estimate','Std. Dev','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
axis tight
xlim([0 27]);
set(gca, 'XTickLabel', []);


%%
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, .4, 1]);
labels_parameters = {'$K_{\mathrm{sat}}$ [m/s]', '$n$ [-]','$\alpha$ [1/m]', '$C_r$ [-]', '$k$ [1/sec]', '$\bar{S}_y$y [-]', '$S_y(1m) [-]$'};
% calibrated_parameters = [7.9433*10^(-5), 1.76, 0.385,-2.33, 0.9, 5e-5];
calibrated_parameters = [1e-6, 1.76, 0.385,-2.33, 0.9, 5.5e-5];
for ii = 1:6
    subplot(8,1,ii + 2)
    if ii < 3
        data_plot = parameterStore(:,:,ii);
    elseif ii <= 5
        data_plot = parameterStore(:,:,ii+1);
    elseif ii <= 6
        data_plot = (Sy_Store_PF)';
    end
    mean_data = mean(data_plot');
    sdt_data = std(data_plot');
    plot(time_plot_days, mean_data, '--','linewidth',0.5,'Color',pallete.blue_colors(2,:));
    % xlim([7,8])
    hold on
    set(gca,'FontSize',12); set(gcf,'Color','White');
    fill([time_plot_days, fliplr(time_plot_days)], [mean_data + sdt_data, fliplr(mean_data - sdt_data)], pallete.blue_colors(3,:), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
    ylabel(labels_parameters{ii},'Interpreter','latex','FontSize',16)
    hold on
    if ii < 3
        plot(time_plot_days,calibrated_parameters(ii)*ones(size(time_plot_days,2)),'--','Color','black','LineWidth',1.5);
    elseif ii < 6
        plot(time_plot_days,calibrated_parameters(ii+1)*ones(size(time_plot_days,2)),'--','Color','black','LineWidth',1.5);
    end     
    % Customize ticks
    set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
    if ii == 1
        % legend('Mean','Std.','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
    end
    xlim([0 27]); 
    if ii ~= 6
        set(gca, 'XTickLabel', []);
    else
        xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); 
    end
end
%saveas(gcf,'Parameters_Only_1.svg')
exportgraphics(gcf,'Parameter_Ensemble.png','Resolution',300)
exportgraphics(gcf,'Parameter_Ensemble.svg', 'ContentType', 'vector');

%% Calculating Sy for a depth
depth_sy = 1; % m
dtheta = parameterStore(:,:,3);
alpha = parameterStore(:,:,4);
n_vg = parameterStore(:,:,2);
Sy_depth = dtheta.*(1 - (1 + (abs(alpha).*depth_sy).^n_vg).^(-(n_vg + 1)./n_vg) ) ;
%% GW Depth
% Load data]
dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\GitHub\2D_Boussinesq_Model\Input_Data_GW_Model_PERTH.xlsx';
input_data = readtable(dir,'Sheet','Assimilation_Data'); % Reading Input data from Excel
z_estimated = table2array(input_data(:,8:9)); % Time and Discharge
% Here you need to specify the observed data
obs_time = round(z_estimated(:,1)); % Time in minutes
discharge = z_estimated(:,2); % Discharge in mm/h

% Gauges Used in the PF
gauges_used = logical(table2array(input_data(20:34,2)));

% Now we go for the GW depth
for j = 1:15 % Up to 15 points
    GW_Depth_Obs(:,j) = table2array(input_data(:,10 + (j-1)*3));
end

% Taking away wells not used
GW_Depth_Obs(:,~gauges_used) = [];

% Now we go for the modeled GW Depths

labels = {'LEO-C_1_0_0_CS451',	'LEO-C_3_0_0_CS451',	'LEO-C_3_1.5_0_CS451',	'LEO-C_3_-1.5_0_CS451',	'LEO-C_3_3.5_0_CS451',	'LEO-C_3_-3.5_0_CS451',	'LEO-C_7_0_0_CS451',	'LEO-C_7_1.5_0_CS451',	'LEO-C_7_-1.5_0_CS451',	'LEO-C_13_0_0_CS451',	'LEO-C_13_1.5_0_CS451',	'LEO-C_13_-1.5_0_CS451',	'LEO-C_17_0_0_CS451',	'LEO-C_21_1.5_0_CS451',	'LEO-C_21_-1.5_0_CS451'};
labels(~gauges_used) = [];
% % Filters
% idx1 = find(time == date_begin);
% idx2 = find(time == date_end);
% time = time(idx1:idx2,:);
% pressure = GW_Depth_Store_PF;
% pressure = pressure(:,:,idx1:idx2);
% pressure_obs = GW_Depth_Obs(idx1:idx2,:);

% % Resampling
% idx_resample = 1:2:length(time);
% time = time(idx_resample);
% pressure = pressure(idx_resample,:); % Now it is every 30 min

% Basic statistics
GW_Depth_Obs(GW_Depth_Obs > 1) = nan; % These are values with wrong simulations stored
GW_Depth_Obs(GW_Depth_Obs < 0) = nan; % These are values with wrong simulations stored
pressure = GW_Depth_Obs(1:size(GW_Depth_Store_PF,3),:);
mean_pressure = mean(pressure);
median_pressure = median(pressure);
std_pressure = std(pressure);
min_pressure = min(pressure);
max_pressure = max(pressure);

% Outliers
GW_Depth_Store_PF(GW_Depth_Store_PF > 1) = nan; % These are values with wrong simulations stored
GW_Depth_Store_PF(GW_Depth_Store_PF < 0) = nan; % These are values with wrong simulations stored

% Plot results
figure('Units', 'normalized', 'Position', [0 0 1 1]);
z_estimated = squeeze(nanmean(GW_Depth_Store_PF,1));
std_estimated_z = squeeze(nanstd(GW_Depth_Store_PF,1));
clear annotation
for i = 1:size(pressure,2)
    subplot(3,ceil(size(pressure,2)/3),i);
    % OBS GW Depth
    plot(time_plot_days, pressure(:,i), 'Color', pallete.red_colors(1,:)); hold on;
    % Modeled GW Depth
    plot(time_plot_days, z_estimated(i,:), '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
    hold on
    set(gca,'FontSize',12); set(gcf,'Color','White');

    % Replace these with your actual data
    time_spinup = 18; % days
    idx = find(time_plot_days >= time_spinup,1,'first');
    zzz1 = pressure(:,i); zzz1(1:idx) = [];
    zzz2 = z_estimated(i,:); zzz2(1:idx) = [];
    [metrics] = calculate_metrics(zzz1, zzz2);

% Add metrics as text annotations on the plot
annotationText = sprintf(['NSE: %.2f\n', ...
                          'RMSE: %.2f\n', ...
                          'PBIAS: %.2f\n', ...
                          'r2: %.2f\n', ...
                          'KGE: %.2f'], ...
                          metrics.NSE, ...
                          metrics.RMSE, ...
                          metrics.PBIAS/100, ...
                          metrics.R_squared, ...
                          metrics.KGE);

% Get current axis position
ax = gca;
pos = ax.Position; % [x, y, width, height]

% Calculate annotation position (normalized units)
x_annot = pos(1) + pos(3) - 0.05; % Slightly inside the right edge
y_annot = pos(2) - 0.1; % Slightly inside the top edge

% Add text box with metrics
annotation('textbox', [x_annot, y_annot, 0.3, 0.25], ...
           'String', annotationText, ...
           'Interpreter', 'latex', ...
           'FontSize', 10, ...
           'EdgeColor', 'none', ...
           'BackgroundColor', 'none','FontName','Montserrat');

    fill([time_plot_days, fliplr(time_plot_days)], [z_estimated(i,:) + std_estimated_z(i,:), fliplr(z_estimated(i,:) - std_estimated_z(i,:))], pallete.blue_colors(3,:), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
    title(labels{i}, 'Interpreter', 'latex');
    ylim([0,1]);
    xlabel('Elapsed Time [days]', 'Interpreter', 'latex'); ylabel('Pressure [m]', 'Interpreter', 'latex');
    if i == 1
        legend({'Measured', 'Particle Filter'}, 'Interpreter', 'latex');
    end

end


%% Pearson Correlation Coefficients between variables
% Load your data
close all
[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS] = coloramps();
% Assuming data has five variables (X1 to X5) and one output (Y - observed discharge)
% Replace these with your actual data
time_spinup = 18; % days
idx = find(time_plot_days >= time_spinup,1,'first');
time = time_plot_days(idx:end); % Example time points (replace with your time series)
zzz = [parameterStore]; 

% Taking away dtheta
zzz(:,:,3) = [];
% Including Sy
zzz(:,:,end+1) = [squeeze(Sy_Store_PF)]';
X = squeeze(mean(zzz(idx:end,:,:),2));
Y = discharge_data(idx:end,1);  % Replace with your observed discharge (output)
Z = GW_Depth_Obs(idx:end-1,:);
Y_PF_mean = nanmean(dischargeStore_PF,2);
Y_PF_mean = Y_PF_mean(idx:end);

% Pearson Correlation
disp('Correlation Matrix:');
corrMatrix = corr([X, Y, Z]); % Include all variables and the output
corrMatrix(eye(size(corrMatrix)) == 1) = NaN;
disp(corrMatrix);

% Create figure
figure('Units', 'centimeters', 'Position', [5, 5, 16, 12]); % Adjust figure size
labels_parameters = {'$K_{\mathrm{sat}}$', '$n$','$\alpha$', '$C_r$', '$k$', '$S_y$'};
labels  = [labels_parameters,'Q', '$W_1$','$W_2$','$W_3$','$W_4$'];
% Display heatmap
h = heatmap(corrMatrix, ...
    'Colormap', cool, ... % Use a nicer colormap (e.g., parula or others)
    'ColorLimits', [-1 1], ... % Correlation limits
    'CellLabelFormat', '%.2f', ... % Format cell values to 2 decimals
    'XDisplayLabels', labels, 'interpreter','latex',...
    'YDisplayLabels', labels, 'interpreter','latex');

% Improve heatmap appearance
h.GridVisible = 'off'; % Remove grid for a cleaner look
h.MissingDataLabel = 'N/A'; % Optional: Handle missing data
h.FontSize = 12; % Increase font size
h.FontName = 'Montserrat'; % Use Montserrat font
h.ColorbarVisible = 'on'; % Add a colorbar

% Adjust figure aesthetics
ax = gca; % Get current axes
% Improve heatmap appearance
h.GridVisible = 'off'; % Remove grid for a cleaner look
h.MissingDataColor = [1 1 1]; % Set missing data (diagonal) to white
h.MissingDataLabel = ''; % Remove label for missing data
h.FontSize = 12; % Increase font size
h.FontName = 'Montserrat'; % Use Montserrat font
h.ColorbarVisible = 'on'; % Add a colorbar
title('Pearson Correlation Coefficients ($r$)')

% exportgraphics(gcf,'Output/Figures/PF/Pearson_Correlation_Matrix.pdf','ContentType','vector');

%%  Spearman Correlation
% Load your data
close all
[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS] = coloramps();
% Assuming data has five variables (X1 to X5) and one output (Y - observed discharge)
% Replace these with your actual data
time_spinup = 18; % days
idx = find(time_plot_days >= time_spinup,1,'first');
time = time_plot_days(idx:end); % Example time points (replace with your time series)
% zzz = [parameterStore]; 
% 
% % Taking away dtheta
% zzz(:,:,3) = [];
% % Including Sy
% zzz(:,:,end+1) = [squeeze(mean(Sy_Store_PF,2))]';
X = squeeze(mean(zzz(idx:end,:,:),2));
Y = discharge_data(idx:end,1);  % Replace with your observed discharge (output)
Z = GW_Depth_Obs(idx:end-1,:);
Y_PF_mean = nanmean(dischargeStore_PF,2);
Y_PF_mean = Y_PF_mean(idx:end);

% 1. Correlation analysis between variables
disp('Correlation Matrix:');
corrMatrix = corr([X, Y, Z],'type','spearman'); % Include all variables and the output
corrMatrix(eye(size(corrMatrix)) == 1) = NaN;
disp(corrMatrix);

% Create figure
figure('Units', 'centimeters', 'Position', [5, 5, 16, 12]); % Adjust figure size
labels_parameters = {'$K_{\mathrm{sat}}$ [m/s]', '$n$ [-]','$\alpha$ [1/m]', '$C_r$ [-]', '$k$ [1/sec]', '$S_y$ [-]'};
labels  = [labels_parameters,'Q', 'W1','W2','W3','W4'];
% Display heatmap
h = heatmap(corrMatrix, ...
    'Colormap', cool, ... % Use a nicer colormap (e.g., parula or others)
    'ColorLimits', [-1 1], ... % Correlation limits
    'CellLabelFormat', '%.2f', ... % Format cell values to 2 decimals
    'XDisplayLabels', labels, 'interpreter','latex',...
    'YDisplayLabels', labels, 'interpreter','latex');

% Improve heatmap appearance
h.GridVisible = 'off'; % Remove grid for a cleaner look
h.MissingDataLabel = 'N/A'; % Optional: Handle missing data
h.FontSize = 12; % Increase font size
h.FontName = 'Montserrat'; % Use Montserrat font
h.ColorbarVisible = 'on'; % Add a colorbar

% Adjust figure aesthetics
ax = gca; % Get current axes
% Improve heatmap appearance
h.GridVisible = 'off'; % Remove grid for a cleaner look
h.MissingDataColor = [1 1 1]; % Set missing data (diagonal) to white
h.MissingDataLabel = ''; % Remove label for missing data
h.FontSize = 12; % Increase font size
h.FontName = 'Montserrat'; % Use Montserrat font
h.ColorbarVisible = 'on'; % Add a colorbar
title('Spearman Correlation Coefficients ($\rho$)')

% exportgraphics(gcf,'Output/Figures/PF/Spearman_Correlation_Matrix.pdf','ContentType','vector');
%% Discharge Comparison
[metrics] = calculate_metrics(Y, Y_PF_mean);
[metrics_calibrated] = calculate_metrics(Y, discharge_benchmark(idx:end));

% Plot observed and simulated time series
figure('Units', 'centimeters', 'Position', [5, 5, 18, 12]); % Figure size
hold on;
plot(time, Y, 'k', 'LineWidth', 2.5, 'DisplayName', 'Observed'); % Observed
hold on
plot(time, Y_PF_mean, 'Color',pallete.red_colors(1,:),'LineStyle','-', 'LineWidth', 2.5, 'DisplayName', 'PF Mean'); % Simulated
hold on
plot(time, discharge_benchmark(idx:end), 'Color',pallete.green_colors(3,:),'LineStyle','--', 'LineWidth', 2.5, 'DisplayName', 'Calibrated'); % Calibrated

% Add plot labels and legend
xlabel('Elapsed Time [days]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Discharge [mm/h]', 'FontSize', 14, 'Interpreter', 'latex');
% title('Hydrograph', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
set(gca,'FontName','Garamond','FontSize',16)
grid on;

% Add metrics as text annotations on the plot
annotationText = sprintf(['NSE: %.2f\n', ...
                          'RMSE: %.2f\n', ...
                          'PBIAS: %.2f%%\n', ...
                          'R^2: %.2f\n', ...
                          'KGE: %.2f'], ...
                          metrics.NSE, ...
                          metrics.RMSE, ...
                          metrics.PBIAS, ...
                          metrics.R_squared, ...
                          metrics.KGE);

% Add text box with metrics
annotation('textbox', [0.15, 0.65, 0.3, 0.25], ...
           'String', annotationText, ...
           'Interpreter', 'latex', ...
           'FontSize', 10, ...
           'EdgeColor', 'none', ...
           'BackgroundColor', 'none','FontName','Montserrat');

% Add metrics as text annotations on the plot
annotationText_calibrated = sprintf(['NSE: %.2f\n', ...
                          'RMSE: %.2f\n', ...
                          'PBIAS: %.2f%%\n', ...
                          'R^2: %.2f\n', ...
                          'KGE: %.2f'], ...
                          metrics_calibrated.NSE, ...
                          metrics_calibrated.RMSE, ...
                          metrics_calibrated.PBIAS, ...
                          metrics_calibrated.R_squared, ...
                          metrics_calibrated.KGE);

% Add text box with metrics
annotation('textbox', [0.35, 0.65, 0.3, 0.25], ...
           'String', annotationText_calibrated, ...
           'Interpreter', 'latex', ...
           'FontSize', 10, ...
           'EdgeColor', 'none', ...
           'BackgroundColor', 'none','FontName','Montserrat');

exportgraphics(gcf,'Output/Figures/PF/Discharge_Comparison.pdf','ContentType','vector');


% Save the figure as a publication-ready image
end

%% Lag-correlation
max_lag = 10; % time-steps
linestyle = {'-','--','-.',':','-','--','-.',':','--','-'};
colors = linspecer((size(X,2) + size(Z,2)));
subplot(1,2,1)
labels_plot = labels; labels_plot(7) = [];
for ii = 1:(size(X,2) + size(Z,2))
if ii <= size(X,2)
    [lags, cross_corr] = correlation_analysis(X(:,ii), Y, max_lag);
else
    [lags, cross_corr] = correlation_analysis(Z(:,ii - size(X,2)), Y, max_lag);
end
plot(lags, cross_corr, 'color',colors(ii,:), 'LineWidth', 2,'LineStyle',linestyle{ii}); % Red line for lag correlation
hold on;
% if ii == 1
%     yline(0, '--k'); % Dashed black zero line
% end
grid on;
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'TickDir', 'out');
xlabel('Lag [30-min intervals]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cross-correlation', 'Interpreter', 'latex', 'FontSize', 14);
title('Lagged Cross-Correlation with Discharge', 'Interpreter', 'latex', 'FontSize', 14);
legend(labels_plot,'interpreter','latex')
hold on;
end
subplot(1,2,2)
for ii = 1:(size(X,2) + 1 + size(Z,2) - 1)
Z_W = Z; Z_W(:,1) = []; % taking away first well
labels_plot = labels; labels_plot(8) = [];
if ii <= size(X,2)
    [lags, cross_corr] = correlation_analysis(X(:,ii), Z(:,1), max_lag);
elseif ii == size(X,2) + 1
    [lags, cross_corr] = correlation_analysis(Y, Z(:,1), max_lag);
else
    [lags, cross_corr] = correlation_analysis(Z_W(:,ii - size(X,2) - 1), Z(:,1), max_lag);
end

plot(lags, cross_corr, 'color',colors(ii,:), 'LineWidth', 2,'linestyle',linestyle{ii}); % Red line for lag correlation
hold on;
% if ii == 1
%     yline(0, '--k'); % Dashed black zero line
% end
end
grid on;
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'TickDir', 'out');
xlabel('Lag [30-min intervals]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cross-correlation', 'Interpreter', 'latex', 'FontSize', 14);
title('Lagged Cross-Correlation with Well 1', 'Interpreter', 'latex', 'FontSize', 14);
legend(labels_plot,'interpreter','latex')
hold on;


%% Error Timeseries
% Plot results
figure('Units', 'normalized', 'Position', [0 0 1 1]);
obs_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Observed_Hydrograph.csv';
z_estimated = table2array(readtable(obs_dir));
z_estimated = z_estimated(:,2);
discharge_data = discharge_data(1:size(dischargeStore_PF,1));
discharge_mean = nanmean(dischargeStore_PF,2);
discharge_std = nanstd(dischargeStore_PF');


z_estimated_PF_mean = squeeze(nanmean(GW_Depth_Store_PF,1));
std_estimated_PF = squeeze(nanstd(GW_Depth_Store_PF,1));
for j = 1:5
    subplot(3,ceil((size(GW_Depth_Store_PF,2) + 1)/3),j);
    if j == 1
    % Discharge
    plot(time_plot_days, abs(discharge_mean - discharge_data), 'Color', pallete.red_colors(1,:)); hold on;
    plot(time_plot_days, sqrt(discharge_std.^2 + 0.1^2), 'Color', pallete.blue_colors(1,:));
    legend('$|\bar{q} - q_{obs}|$','$\sqrt{\sigma_{\mathbf{q}}^2 + \sigma_{obs}^2}$','interpreter','latex')
    xlabel('Elapsed time [day]','interpreter','latex')
    ylabel('Error ','interpreter','latex')
    else
    i = j -1;
    % Obs Depth
    GW_depth_std = squeeze(GW_Depth_Store_PF(:,i,:));
    GW_depth_std = nanstd(GW_depth_std);
    plot(time_plot_days, abs(z_estimated_PF_mean(i,:) - GW_Depth_Obs(1:1296,i)'), 'Color', pallete.red_colors(1,:)); hold on;
    plot(time_plot_days, sqrt(GW_depth_std.^2 + 0.05^2)', 'Color', pallete.blue_colors(1,:)); hold on;
    xlabel('Elapsed time [day]','interpreter','latex')
    ylabel('Error ','interpreter','latex')
    end
end


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
