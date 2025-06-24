%% Clear and Load Data
close all; clear; clc;

clear all
load 'Output/PERTH_CALIBRATION_400_particles_smaller';

% Load observed hydrograph
obs_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Observed_Hydrograph.csv';
obs_data = table2array(readtable(obs_dir));
discharge_data = obs_data(:,2); % Observed discharge

% Load calibrated (benchmark) hydrograph
benchmark_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Calibration\calibrated_discharge.csv';
cal_data = table2array(readtable(benchmark_dir));
discharge_benchmark = cal_data(:,2); % Calibrated discharge

% Ensure size consistency
discharge_data = discharge_data(1:size(dischargeStore_PF,1));
discharge_benchmark = discharge_benchmark(1:size(dischargeStore_PF,1));

%% Compute PF Estimated Discharge
dischargeStore_PF(dischargeStore_PF>100) = nan;
z_estimated = nanmean(dischargeStore_PF, 2); % Mean discharge across realizations
std_estimated_z = nanstd(dischargeStore_PF, 0, 2); % Standard deviation

%% Time and Rainfall Setup
dt = 0.5; % Time step in hours
time = (0:length(discharge_data)-1) * dt; % Time vector

% Extract Rainfall (TopBC is 3D: [time, row, column])
rainfall = [11.73577588	11.58073279	11.46716153	11.00498194	10.80300909	10.81407012	3.356556921	0	0	0	0	0	0	0	0	0	0	0	0	0	11.22076621	11.48278151	11.47776421	11.46415795	11.52082143	11.51993467	1.556299634	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	8.756336469	10.28481552	11.44924116	10.88548889	21.24002009	23.24649304	6.251009065	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.66075022	12.03314426	11.80620709	11.59077604	0.823002711	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.19151164	11.74289861	11.65096982	11.33082764	11.53962188	11.71976012	1.910525576	0	0	0	0	0	0	0	0	0	0	0	0	0	11.28669394	11.57983048	11.55855152	11.61208818	11.57156685	11.52564612	1.1015656	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.21932641	11.64994195	11.66286644	11.68306313	11.67958496	11.64166745	0.670342924	0	0	0	0	0	0	0	0	0	0	0	0	0	11.25192315	11.61933176	11.66508873	11.64706721	11.53584236	11.62283291	1.881064348	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.28064442	11.78186109	11.70809055	11.39201527	11.57706461	11.48519582	1.517066923	0	0	0	0	0	0	0	0	0	0	0	0	0	11.38635091	11.73607376	11.6793517	11.64315685	11.65822103	11.68336382	0.645044121	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.42739218	11.86222255	11.73671006	11.74549885	11.71035145	11.59079745	0.161112485	0	0	0	0	0	0	0	0	0	0	0	0	0	11.21472467	11.65388097	11.7070323	11.74669012	11.71311861	11.72362806	0.68350138	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.28754721	11.70993333	11.73509145	11.78277958	11.72465752	11.63938121	0.913347879	0	0	0	0	0	0	0	0	0	0	0	0	0	11.27684618	11.78238194	11.79905012	11.77841212	11.70930861	11.69559115	1.036390758	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11.29139891	11.74286624	11.75372958	11.69981952	11.76092352	11.69026376	0.93455001	0	0	0	0	0	0	0	0	0	0	0	0	0	11.00109782	11.8150543	11.7614617	11.77735158	11.74633703	11.70717655	0.785077564	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
rainfall = rainfall(1,1:end-1)';
% Convert rainfall to m³/s assuming rainfall over 1 m² (1 mm = 1 L/m² = 1e-3 m³/m²)
% If you have catchment area, multiply rainfall by area to get total inflow.
rainfall_m3h = rainfall / 1000 * 60 * 24 * 0.47; % m/s to m³/h over area

%% Spinup array
days = 14;
nt = days / (dt/24);
discharge_data(1:nt) = [];
z_estimated(1:nt) = [];
rainfall(1:nt) = [];
rainfall_m3h(1:nt) = [];
discharge_benchmark(1:nt) = [];
time(1:nt) = [];

%% Compute Outflow Volumes
PFoutVol = nansum(z_estimated * dt); % Simulated total outflow (m³)
obsOutFlow = nansum(discharge_data * dt); % Observed total outflow (m³)
obsCalFlow = nansum(discharge_benchmark * dt); % Calibrated total outflow (m³)

%% Compute PBIAS
PBIAS_obs = 100 * (PFoutVol - obsOutFlow) / obsOutFlow; % PBIAS vs Observed
PBIAS_cal = 100 * (PFoutVol - obsCalFlow) / obsCalFlow; % PBIAS vs Calibrated

% Display PBIAS results
fprintf('PBIAS with Observed: %.2f%%\n', PBIAS_obs);
fprintf('PBIAS with Calibrated: %.2f%%\n', PBIAS_cal);

%% Compute Storage Over Time
% Cumulative rainfall input (m³)
cum_rain_input = cumsum(rainfall_m3h * dt);

% Cumulative outflows (m³)
cum_outflow_PF = cumsum(z_estimated * dt);
cum_outflow_obs = cumsum(discharge_data * dt);
cum_outflow_cal = cumsum(discharge_benchmark * dt);

% Storage over time
storage_PF = cum_rain_input - cum_outflow_PF;
storage_obs = cum_rain_input - cum_outflow_obs;
storage_cal = cum_rain_input - cum_outflow_cal;

%% Plot Hydrographs and Rainfall
figure;

yyaxis left
plot(time, discharge_data, 'b-', 'LineWidth', 1.5); hold on;
plot(time, discharge_benchmark, 'g--', 'LineWidth', 1.5);
plot(time, z_estimated, 'r-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Discharge (m^3/s)');
legend('Observed', 'Calibrated', 'PF Simulated', 'Location', 'best');
grid on;

yyaxis right
bar(time, rainfall, 0.5, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Rainfall as bar plot
ylabel('Rainfall (mm/h)');

title('Hydrograph and Rainfall');
set(gca, 'FontSize', 12);

%% Plot Storage Comparison
figure;
plot(time, storage_obs, 'b-', 'LineWidth', 1.5); hold on;
plot(time, storage_cal, 'g--', 'LineWidth', 1.5);
plot(time, storage_PF, 'r-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Storage (m^3)');
legend('Observed Storage', 'Calibrated Storage', 'PF Simulated Storage', 'Location', 'best');
grid on;
title('Storage Over Time');
set(gca, 'FontSize', 12);
