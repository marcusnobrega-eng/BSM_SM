%% Testing Various Cases
clear all

% Model Parameters
% Ksat,h: horizontal ksat [m/s]
% n_values: UZ Van Genutchen alpha values [1/m]
% Alpha_values: UZ Van Genutchen alpha values [-]
% k: linear reservoir factor for the UZ damping [1/sec] (1e-4 typically)
% r_damping: recharge damping factor Irrigation = Irrigation*r_damping

% Ksat_values = [5e-5:1e-5:5e-4];
% n_values = [1.5:0.1:2.5];
% Poro_values = [0.37];
% Alpha_values = [-2.33];

% Ksat_values = [1e-4:1e-5:4e-4];
Ksat_values = 10.^(-[3:0.1:6]);
n_values = 1.76;
Poro_values = 0.385;
Alpha_values = -2.33;
k_values = [0.5e-5:0.5e-5:1e-4];
irr_damping = [0.90];

% x = [Ksat, n, Poro, Alpha, k]

% Ksat_values = [7.9109e-05];
% n_values = 4.93;
% Poro_values = 0.385;
% Alpha_values = -4.08;
% 
% % Automatic Calibration Resuls
% Ksat_values = 6.9e-05;
% n_values = 4.93;
% Poro_values = 0.37;
% Alpha_values = -4.08;

% n_values = [0.3];
% Ksat_values = [0.75*10^(-4)];

% n_values = [0.35];
% Ksat_values = 7e-5;

dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\GitHub\2D_Boussinesq_Model\Input_Data_GW_Model_PERTH.xlsx';
index = 1;
flag_opt_fun = 1;
n_tests = length(Ksat_values)*length(k_values);
discharge_values = zeros(1296,length(k_values)*length(Ksat_values));
GW_values = zeros(1296,4,n_tests);

% i = 13, j = 7
for j = 1:length(Ksat_values)
    for i = 1:length(k_values)
        % Optimal Values i = 10; j = 12; Ksat = 7.9433e-5 m/s, k = 5e-5
        % 1/sec
        x = [Ksat_values(j), n_values(1), Poro_values(1), Alpha_values(1), k_values(i), irr_damping(1)];
        [Obj_fun,discharge,GW_Depth_Modeled] = Opt_function(x, dir, flag_opt_fun);
        discharge_values(:,index) = discharge';
        GW_values(:,:,index) = GW_Depth_Modeled;
        Performance_values(index,1) = -Obj_fun;
        Performance_matrix(i,j) = -Obj_fun;
        index = index + 1;
    end
end

%% Plotting Results
load 'calibrated_results_GW_discharge'
[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS] = coloramps();
% load calibrated_data_2.mat;
input_table = readtable('Observed_Hydrograph.csv');
Q_obs_table = table2array(input_table(:,:));
time_Q_obs = Q_obs_table(:,1);
time_Q_obs = time_Q_obs/1440; % days
Q_obs = Q_obs_table(:,2:end); % only last entries

figure(3)
% Load irrigation data
input_table = readtable('Calibration\Irrigation_data.csv');
I_obs_table = table2array(input_table(:,:));
time_I_obs = I_obs_table(:,1);
time_I_obs = time_I_obs/1440; % days
I_obs = I_obs_table(:,2:end); % only last entries

% Summing all irrigation
irr_tot = sum(I_obs*(time_I_obs(2) - time_I_obs(1)))*24; % mm

% Summing all outflow
outflow_tot = sum(Q_obs*(time_Q_obs(2) - time_Q_obs(1))*24);
zzz = Performance_values;
zzz(isinf(zzz)) = nan;
idx = find(zzz == max(max(zzz)));
Performance_values(idx)
modoutflow_tot = sum(1000*3600*discharge_values(:,idx)*(time_Q_obs(2) - time_Q_obs(1))*24);
GW_Depth_Optimal = GW_values(:,:,idx);

%% Optimal Values
% Define the index for which you want to retrieve x
target_index = idx;

% Assuming you know the length of the loops
num_k_values = length(k_values);
num_Ksat_values = length(Ksat_values);

% Calculate the corresponding (i, j) values
j = ceil(target_index / num_k_values);  % The outer loop index
i = target_index - (j - 1) * num_k_values;  % The inner loop index

% Reconstruct x based on the values of i and j
x_reconstructed = [Ksat_values(j), n_values(1), Poro_values(1), Alpha_values(1), k_values(i), irr_damping(1)];

%% GW Depth
t_step_begin = 673; % Number of steps to start
input_data = readtable(dir,'Sheet','Assimilation_Data'); % Reading Input data from Excel
z_estimated = table2array(input_data(:,8:9)); % Time and Discharge
% Here you need to specify the observed data

cmf = General_Setup_Calibrator(x,dir);

% Now we go for the GW depth
DX = cmf.coupledModel.gwmodel.nx*cmf.coupledModel.gwmodel.dx;
DY = cmf.coupledModel.gwmodel.ny*cmf.coupledModel.gwmodel.dy;
for j = 1:15 % Up to 15 points
    GW_Depth_Obs(:,j) = table2array(input_data(:,10 + (j-1)*3));
    x_well(j,1) = ceil(table2array(input_data(1,11 + (j-1)*3))/cmf.coupledModel.gwmodel.dx); % Number of cells
    y_well(j,1) = ceil((DY - table2array(input_data(1,12 + (j-1)*3)))/cmf.coupledModel.gwmodel.dy); % Number of cells
end

GW_Depth_Obs(1,:) = []; % First value not used in the calibration

% Gauges Used in the PF
gauges_used = logical(table2array(input_data(20:34,2)));

z_theoretical_depth = GW_Depth_Obs;
% z_theoretical_depth = GW_Depth_Obs(t_step_begin:end,:);
% 
% GW_Depth_Obs = GW_Depth_Obs(t_step_begin:end,:);

% Taking away wells not used
z_theoretical_depth(:,~gauges_used) = [];
x_well(~gauges_used) = [];
y_well(~gauges_used) = [];

Optimal_Discharge = discharge_values(:,idx);
%% Plotting Best Case
plot(time_Q_obs,Q_obs,'LineWidth',2,'color',pallete.blue_colors(2,:));
hold on
plot(time_Q_obs,1000*3600*discharge_values(:,idx),'LineWidth',2,'color',pallete.red_colors(2,:));
xlabel('Elapsed time [days]','Interpreter','latex');
ylabel('Discharge [mm/h]','Interpreter','latex');
ylim([0,2])
yyaxis right
set(gca,'ydir','reverse');
set(gca,'YColor','black')
ylabel('Irrigation [mm/h]','Interpreter','latex');
bar([0;time_Q_obs],I_obs,'FaceColor',[0 128 128]/256,'EdgeColor',[.3 .3 .3],'LineWidth',1.2)
ylim([0,50]);
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
title('Calibrated Model','Interpreter','latex','FontSize',14);
legend('Observed','Best','Irrigation','interpreter','latex','Location','best')
exportgraphics(gcf,'calibration_discharge.pdf','Resolution',600)

%% Plotting Best Case for GW
for i = 1:4
subplot(2,2,i)
plot(time_Q_obs,GW_Depth_Obs(:,i),'LineWidth',2,'color',pallete.blue_colors(2,:));
hold on
plot(time_Q_obs,GW_Depth_Optimal(:,i),'LineWidth',2,'color',pallete.red_colors(2,:));
xlabel('Elapsed time [days]','Interpreter','latex');
ylabel('Discharge [mm/h]','Interpreter','latex');
ylim([0,1])
yyaxis right
set(gca,'ydir','reverse');
set(gca,'YColor','black')
ylabel('Irrigation [mm/h]','Interpreter','latex');
bar([0;time_Q_obs],I_obs,'FaceColor',[0 128 128]/256,'EdgeColor',[.3 .3 .3],'LineWidth',1.2)
ylim([0,50]);
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
title('Calibrated Model','Interpreter','latex','FontSize',14);
if i == 1
legend('Observed','Best','Irrigation','interpreter','latex','Location','best')
end
end
exportgraphics(gcf,'calibration_discharge.pdf','Resolution',600)

%%
plot(Ksat_values(1:length(Performance_values)),Performance_values,'LineWidth',2,'color',pallete.blue_colors(2,:));
xlabel('$K_{sat,h}$ [m/s]','Interpreter','latex');
ylabel('NSE [-]','Interpreter','latex');


%% Performance Surface
% load calibrated_data_dampen.mat
[X,Y] = meshgrid(Ksat_values,k_values);
cut_performance = Performance_matrix;
cut_performance(cut_performance<-1) = nan;
hold on
figure(1)
 [~] = surf_plot_full(max(max(cut_performance)),10,'\mathrm{NSE}','-',cut_performance,0,0,100,0.9,0,[0,90],X,Y,'$K_{\mathrm{sat}}$ [$\mathrm{m \cdot s^{-1}}$]','$k$ [$\mathrm{sec}^{-1}$]');
 xscale log
 yscale log
axis tight
axis equal
exportgraphics(gcf,'calibration.pdf','Resolution',300)
%%
figure(3)
figure;
[contour_plot, h] = contour(X, Y, cut_performance, 30, 'LineWidth', 1.5);  % 15 contour levels
% h.LevelList = round(h.LevelList,2);  % rounds levels to 3rd decimal place
% Add labels to the contour lines with 2 decimal places
% clabel(contour_plot, h, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black', 'LabelSpacing', 500);  % Use two decimal places for contour labels

% Set the font for the entire plot to LaTeX
set(gca, 'FontName', 'Helvetica');  % Helvetica (default for LaTeX)
set(gca, 'TickLabelInterpreter', 'latex');  % Ensure ticks are in LaTeX font
xlabel('$K_{\mathrm{sat,h}}$', 'Interpreter', 'latex','FontSize',12);
ylabel('$k$', 'Interpreter', 'latex', 'Interpreter', 'latex','FontSize',12);

% Add colorbar with turbo colormap
colorbar_handle = colorbar;
colormap turbo;  % Use turbo colormap
ylabel(colorbar_handle, 'Performance', 'Interpreter', 'latex');  % Label the colorbar

% Add ticks on the outside
set(gca, 'TickDir', 'out');  % Ticks on the outside
set(gca, 'Box', 'on');       % Ensure the box is around the plot
xscale log
yscale log
axis tight
% Show the plot
grid on;
exportgraphics(gcf,'decision_space.pdf','Resolution',600)

figure(2)
index = 1;
colors = linspecer(length(Ksat_values));
for i = 1:length(n_values)
    subplot(length(n_values),1,i)
    plot(time_Q_obs,Q_obs,'LineWidth',2);
    hold on
    for j = 1:length(Ksat_values)        
        plot(time_Q_obs,1000*3600*discharge_values(:,index),'LineWidth',1.3,'Color',colors(j,:),'LineStyle','--');
        index = index + 1;
        hold on
    end
    xlabel('Elapsed time [days]');
    ylabel('Discharge [mm/h]');
end

figure(3)
% Load irrigation data
input_table = readtable('Calibration\Irrigation_data.csv');
I_obs_table = table2array(input_table(:,:));
time_I_obs = I_obs_table(:,1);
time_I_obs = time_I_obs/1440; % days
I_obs = I_obs_table(:,2:end); % only last entries

% Summing all irrigation
irr_tot = sum(I_obs*(time_I_obs(2) - time_I_obs(1)))*24; % mm

% Summing all outflow
outflow_tot = sum(Q_obs*(time_Q_obs(2) - time_Q_obs(1))*24);
zzz = Performance_values;
zzz(isinf(zzz)) = nan;
idx = find(zzz == max(max(zzz)));
modoutflow_tot = sum(1000*3600*discharge_values(:,idx)*(time_Q_obs(2) - time_Q_obs(1))*24);

% Plotting Best Case
plot(time_Q_obs,Q_obs,'LineWidth',2,'color',pallete.blue_colors(2,:));
hold on
plot(time_Q_obs,1000*3600*discharge_values(:,idx),'LineWidth',2,'color',pallete.red_colors(2,:));
xlabel('Elapsed time [days]','Interpreter','latex');
ylabel('Discharge [mm/h]','Interpreter','latex');
ylim([0,2])
yyaxis right
set(gca,'ydir','reverse');
set(gca,'YColor','black')
bar([0;time_Q_obs],I_obs,'FaceColor',[0 128 128]/256,'EdgeColor',[.3 .3 .3],'LineWidth',1.2)
ylim([0,50]);
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
title('Calibrated Model','Interpreter','latex','FontSize',14);