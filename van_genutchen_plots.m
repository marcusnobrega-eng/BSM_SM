% Van Genuchten Parameters for Soil (example values)
calibrated_parameters = [7.9433*10^(-5), 1.76, 0.385,-2.33, 0.9, 5e-5];
theta_r = 0.0;  % Residual water content (cm^3/cm^3)
theta_s = 0.385;  % Saturated water content (cm^3/cm^3)
alpha = 2.33;    % Alpha parameter (1/cm)
n = 1.75;         % n parameter (dimensionless)
m = 1 - 1/n;     % m parameter (dimensionless)
n_points = 100;

[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps();

% Range of water content (x-axis) from residual to saturated
theta = linspace(theta_r + 0.01, theta_s, n_points);  % Water content (cm^3/cm^3)

% Van Genuchten Equation for Soil Suction Head (h)
h = ( (theta_s - theta_r) ./ (theta - theta_r) ).^(1/n) - 1;
h = h .^ (1/m) / alpha;  % Soil suction head in cm

% Particle Filter Results (example n and alpha values)
time = 24.5; % days
Nt = round(time/0.5); % time-steps
n_values = squeeze(parameterStore(Nt,:,2));  % Example particle filter results for n (e.g., 10 values)
alpha_values = squeeze(abs(parameterStore(Nt,:,4)));  % Example particle filter results for alpha (e.g., 10 values)

% Initialize arrays to store the ensemble results
mean_alpha = nanmean(alpha_values);
mean_n = nanmean(n_values);
std_alpha = nanstd(alpha);
std_n = nanstd(n);

h_PF = zeros(length(n_values),n_points);

% Calculate the Van Genuchten equation for each particle filter result
for i = 1:length(n_values)
    n_PF = n_values(i);
    alpha_PF = alpha_values(i);
    
    % Van Genuchten Equation for Soil Suction Head (h) for each n and alpha
    h_PF(i,:) = ( (theta_s - theta_r) ./ (theta - theta_r) ).^(1/n_PF) - 1;
    h_PF(i,:) = h_PF(i,:) .^ (1/(1 - 1/n_PF)) / alpha_PF;  % Soil suction head in cm
   
end
h_PF(isinf(h_PF)) = nan;
% Add to ensemble mean and standard deviation
ensemble_mean = nanmean(h_PF);
ensemble_std = nanstd(h_PF);

% Upper and lower bounds of the ensemble
ensemble_upper = ensemble_mean + ensemble_std;
ensemble_lower = ensemble_mean - ensemble_std;

% Plotting
figure;
plot(theta, h, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1]);  % Orange color
set(gca, 'YScale', 'log');  % Log scale for y-axis
set(gca, 'TickDir', 'out');  % Ticks outside
hold on;

% Plot the ensemble mean
plot(theta, ensemble_mean, 'LineWidth', 2, 'Color', [0, 0.5, 0]);  % Green color for mean
hold on

plot(theta, ensemble_upper, 'LineWidth', 2, 'LineStyle','-.','Color', pallete.blue_colors(3,:));

hold on

plot(theta, ensemble_lower, 'LineWidth', 1.5, 'LineStyle','-.','Color', pallete.blue_colors(3,:));  % 

legend('Calibrated','PF Mean','Upper Boundary','Lower Boundary','interpreter','latex','FontSize',14);

% Plot the ensemble upper and lower bounds as shaded area
% fill([theta, fliplr(theta)], [ensemble_lower, fliplr(ensemble_upper)], ...
%     'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % Green fill for bounds with transparency

% fill([theta, fliplr(theta)], [ensemble_upper, fliplr(ensemble_lower)], pallete.blue_colors(3,:), ...
%         'FaceAlpha', 0.3, 'EdgeColor', 'none');

% set(gca, 'YScale', 'log');  % Log scale for y-axis

% Fill the area with 'patch' for log-scale compatibility
% hold on
% X = [theta, fliplr(theta)];
% Y = [ensemble_upper, fliplr(ensemble_lower)];
% patch(X, Y, [0.7, 1, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Light green shading

% % Plot ensemble mean
% semilogy(theta, ensemble_mean, 'g', 'LineWidth', 2);  % Green line
% 
% % Plot Van Genuchten curve using the first n and alpha values
% semilogy(theta, h(1, :), 'r--', 'LineWidth', 2);  % Red dashed line for reference


% Labeling with LaTeX formatting
xlabel('$\theta$ $[\mathrm{cm^3/cm^3}]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Suction Head ($\psi$) [cm]', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('t = %.2f days', time), ...
    'Interpreter', 'latex', 'FontSize', 16);
% Customizing grid and figure properties
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.7, 0.7, 0.7]);  % Light gray grid lines

% Set tick font size and style
set(gca, 'FontSize', 14, 'FontName', 'Garamond', 'TickLength', [0.02, 0.03]);

% Optional: You can change axis limits if needed for better view
xlim([theta_r, theta_s]);

% Now add particle filter results

