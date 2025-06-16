% Pressure Analysis
close all
[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS] = coloramps();

% Dates
date_begin = datetime(2016,12,1,8,0,0);
date_end = datetime(2016,12,28,8,0,0);

% Load data
data = readtable('pressure_data_PERTH.xlsx'); % Assumes first column is time, next 15 are pressure sensors
time = data{3:end,1};
pressure = data{3:end,2:end};
pressure(pressure<0) = 0;

labels = {'LEO-C_1_0_0_CS451',	'LEO-C_3_0_0_CS451',	'LEO-C_3_1.5_0_CS451',	'LEO-C_3_-1.5_0_CS451',	'LEO-C_3_3.5_0_CS451',	'LEO-C_3_-3.5_0_CS451',	'LEO-C_7_0_0_CS451',	'LEO-C_7_1.5_0_CS451',	'LEO-C_7_-1.5_0_CS451',	'LEO-C_13_0_0_CS451',	'LEO-C_13_1.5_0_CS451',	'LEO-C_13_-1.5_0_CS451',	'LEO-C_17_0_0_CS451',	'LEO-C_21_1.5_0_CS451',	'LEO-C_21_-1.5_0_CS451'};

% Filters
idx1 = find(time == date_begin);
idx2 = find(time == date_end);
time = time(idx1:idx2,:);
pressure = pressure(idx1:idx2,:);

% Resampling
idx_resample = 1:2:length(time);
time = time(idx_resample);
pressure = pressure(idx_resample,:); % Now it is every 30 min

% Basic statistics
mean_pressure = mean(pressure);
median_pressure = median(pressure);
std_pressure = std(pressure);
min_pressure = min(pressure);
max_pressure = max(pressure);

% Noise detection using moving standard deviation
window_size = 10; % Adjust for sensitivity
moving_std = movstd(pressure, window_size);
thresh = mean(moving_std) + 4*std(moving_std); % Define noise threshold
noisy_points = moving_std > thresh;

% Smoothness check using moving average
smoothed_pressure = movmean(pressure, window_size);

% Set LaTeX font properties
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultAxesFontSize', 15);
set(groot, 'DefaultAxesTickDir', 'out');

% Define colors
colors = lines(15);

% Plot results
figure('Units', 'normalized', 'Position', [0 0 1 1]);

for i = 1:15
    subplot(3,5,i);
    plot(time, pressure(:,i), 'Color', pallete.red_colors(1,:)); hold on;
    scatter(time(noisy_points(:,i)), pressure(noisy_points(:,i),i), 'Color','k')
    title(labels{i}, 'Interpreter', 'latex');
    xlabel('Date', 'Interpreter', 'latex'); ylabel('Pressure [m]', 'Interpreter', 'latex');
    if i == 1
    legend({'Raw Data', 'Noisy Points'}, 'Interpreter', 'latex');
    end
end

figure('Units', 'normalized', 'Position', [0 0 1 1]);

for i = 1:15
    subplot(3,5,i);
    plot(time, moving_std(:,i), 'k'); hold on;
    yline(thresh(i), 'Color',pallete.red_colors(1,:),'LineStyle','--');
    title(labels{i}, 'Interpreter', 'latex');
    xlabel('Date', 'Interpreter', 'latex'); ylabel('Std Dev [m]', 'Interpreter', 'latex');
    if i == 1
    legend({'Moving Std', 'Threshold'}, 'Interpreter', 'latex');
    end
end

figure('Units', 'normalized', 'Position', [0 0 1 1]);

for i = 1:15
    subplot(3,5,i);
    plot(time, pressure(:,i), 'Color',pallete.blue_colors(1,:)); hold on;
    plot(time, smoothed_pressure(:,i), 'Color',pallete.red_colors(1,:));
    title(labels{i}, 'Interpreter', 'latex');
    xlabel('Date', 'Interpreter', 'latex'); ylabel('Pressure [m]', 'Interpreter', 'latex');
    if i == 1
    legend({'Raw Data', 'Smoothed Data'}, 'Interpreter', 'latex');
    end
end

disp('Analysis complete. Check the plots for noise detection and smoothness evaluation.');
