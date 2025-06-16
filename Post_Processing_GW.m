% Post-Processing Code for GW Model
% Developer: Marcus Nobrega, Ph.D.
% Last Updated: 4/5/2024

%% Load Data
% load 'Output'/LEO_test.mat
% load 'Output'/Bucket.mat
% load 'Output'/miniLEO.mat
% load 'Output'/LEO.mat
% load 'GW_10perc_slope_drainage_20perc'
% load 'Output'/GW_5_perc_drainage;
% load 'Output'/GW_5_perc_recharge;
% load 'Output'/GW_30_perc_drainage;
% load 'Output'/GW_30_perc_drainage;
% load 'Output'/GW_30_perc_recharge;
% load 'Output'/GW_30_perc_drainage;
% load 'Output'/GW_5_perc_drainage;
% load GW_10_perc_recharge_convergent
% load GW_10_perc_drainage_convergent
% load GW_10_perc_recharge_convergent

close all

%% Addpath where plotting functions are located
path_fun = 'C:\Users\Marcus\Desktop\Desktop_Folder\LEO\h3d-brandhorst-erdal-main-CoupledIterative\CoupledIterative\Extra';
try
    addpath path_fun
catch 
    error('Please add the path where the output functions are located.')
end
close all

%% Creating Folders
try 
    mkdir Output
catch 
    warning('Folder (Output) already created')
end

try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end

try 
    mkdir Output/Figures
catch 
    warning('Folder (Figures) already created')
end

%% Coloramp
[Spectrum,depth_ramp,terrain_ramp] = coloramps(); % Run coloramp function

%% 2.0 Data
time_step = tSave(1,1);                 % Seconds
time_plot_days = tSave/(3600*24);       % days

%% 3.0 Reshaping Data to 2D matrices
nx = size(lsurf,2); % [m]
ny = size(lsurf,1); % [m]
[X,Y] = meshgrid(dx:dx:nx*dx,dy:dy:dy*ny);

% Groundwater Levels
GW_table = GW + z0; % m

for i = 1:size(GW,3)
    matrix = GW(:,:,i);
    GW_Storage(i,1) = sum(sum(matrix(~isnan(matrix))))*dx*dy; % m3
end

%% Colors
n_plots = 10;
colors = linspecer(n_plots,1);

% Outlet Cell
[col_out,row_out] = find(lsurf == min(min(lsurf)),1,'first');

% Internal Cell
[col_int,row_int] = find(lsurf >= max(max(lsurf))/2,1,'first');

%% Bedrock Elevation
% elevation_surf(time_plot_days,'Bedrock Elevation','z','m',z0,X,Y)

%% Only discharge
flag_PF = 1;
if flag_PF == 1
% Particle Filter Data
dischargeStore_PF(dischargeStore_PF > 100) = nan; % These are values with wrong simulations stored

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
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
plot(time_plot_days, discharge_data, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_z, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(time_plot_days, z_estimated, '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([time_plot_days, fliplr(time_plot_days)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); ylabel('Discharge [mm/h]','Interpreter','latex','FontSize',16)
plot(time_plot_days, discharge_benchmark, '-','linewidth',1.5,'Color',pallete.red_colors(2,:));
% xlim([7,8])
ylim([0,2]);
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Observed','Particle filter estimate','Std. Dev','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
axis tight
xlim([0 27]);
exportgraphics(gcf,'Discharge_Ensemble.png','Resolution',300)
end

if flag_PF == 1
    error('Stop here otherwise it will take too long')
end
%% Particle Filter Plots
flag_PF = 1;
if flag_PF == 1
% Particle Filter Data
dischargeStore_PF(dischargeStore_PF > 100) = nan; % These are values with wrong simulations stored

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
subplot(4,2,[1])
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
plot(time_plot_days, discharge_data, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_z, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(time_plot_days, z_estimated, '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([time_plot_days, fliplr(time_plot_days)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); ylabel('Discharge [mm/h]','Interpreter','latex','FontSize',16)
plot(time_plot_days, discharge_benchmark, '-','linewidth',1.5,'Color',pallete.red_colors(2,:));
% xlim([7,8])
ylim([0,2]);
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Observed','Particle filter estimate','Std. Dev','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
axis tight
xlim([0 27]);

subplot(4,2,[2])
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
plot(time_plot_days, discharge_data, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_z, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(time_plot_days, z_estimated, '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([time_plot_days, fliplr(time_plot_days)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); ylabel('Discharge [mm/h]','Interpreter','latex','FontSize',16)
plot(time_plot_days, discharge_benchmark, '-','linewidth',1.5,'Color',pallete.red_colors(2,:));
% xlim([7,8])
ylim([0,2]);
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Observed','Particle filter estimate','Std. Dev','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
axis tight
xlim([0 27]);

hold on
labels_parameters = {'$K_{sat,h}$ [m/s]', '$n$ [-]','$\Delta \theta$ [-]','$\alpha$ [1/m]', '$i_r$ [-]', '$k$ [1/sec]'};
calibrated_parameters = [7.9433*10^(-5), 1.76, 0.385,-2.33, 0.9, 5e-5];
for ii = 1:6
    subplot(4,2,ii + 2)
    data_plot = parameterStore(:,:,ii);
    mean_data = mean(data_plot');
    sdt_data = std(data_plot');
    plot(time_plot_days, mean_data, '--','linewidth',0.5,'Color',pallete.blue_colors(2,:));
    % xlim([7,8])
    hold on
    set(gca,'FontSize',12); set(gcf,'Color','White');
    fill([time_plot_days, fliplr(time_plot_days)], [mean_data + sdt_data, fliplr(mean_data - sdt_data)], pallete.blue_colors(3,:), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
    xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); 
    ylabel(labels_parameters{ii},'Interpreter','latex','FontSize',16)
    hold on
    plot(time_plot_days,calibrated_parameters(ii)*ones(size(time_plot_days,2)),'--','Color','black','LineWidth',1.5);
    % if ii == 1
    %     yscale log
    % end
    % Customize ticks
    set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
    legend('$\mu$','$\sigma$','Calibrated','Interpreter','latex','FontSize',16,'Location','Best')
    xlim([0 27]);    
end
exportgraphics(gcf,'Discharge_Ensemble_7_8.png','Resolution',300)
end

%% Particle Filter Plots
flag_PF = 1;
if flag_PF == 1
close all
% load z_estimated.csv
obs_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Observed_Hydrograph.csv';
z_estimated = table2array(readtable(obs_dir));

obs_time = round(z_estimated(:,1));
discharge_data = z_estimated(:,2);
% discharge_data = discharge_data(1:size(dischargeStore_PF,1));

bechmark_dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Calibration\calibrated_discharge.csv';
z_estimated = table2array(readtable(obs_dir));
discharge_benchmark = z_estimated(:,2); % benchmark discharge

z_estimated = mean(dischargeStore_PF');
z_median = median(dischargeStore_PF');
z_estimated = z_median;
std_estimated_z = std(dischargeStore_PF');

[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps();
figure(1);
subplot(4,1,1)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
plot(obs_time/1440, discharge_data, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_z, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(time_plot_days, z_estimated, '--','linewidth',.5,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([time_plot_days, fliplr(time_plot_days)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); ylabel('Discharge [mm/h]','Interpreter','latex','FontSize',16)
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Observed','Particle filter estimate','Interpreter','latex','FontSize',16,'Location','Best')
ylim([0,3])
hold on
labels_parameters = {'$K_{sat,h}$ [m/s]','$S_y$ [-]'};
for ii = 1:2
    subplot(4,1,ii + 1)
    data_plot = squeeze(parameterStore(:,:,ii));
    mean_data = mean(data_plot');
    sdt_data = std(data_plot');
    plot(time_plot_days, mean_data, '--','linewidth',0.5,'Color',pallete.blue_colors(2,:));
    hold on
    set(gca,'FontSize',12); set(gcf,'Color','White');
    fill([time_plot_days, fliplr(time_plot_days)], [mean_data + sdt_data, fliplr(mean_data - sdt_data)], pallete.blue_colors(3,:), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
    xlabel('Elapsed time [days]','Interpreter','latex','FontSize',16); 
    ylabel(labels_parameters{ii},'Interpreter','latex','FontSize',16)
    % if ii == 1
    %     set(gca,'YScale','log');
    % end
    % Customize ticks
    set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
    legend('$\mu$','$\sigma$','Interpreter','latex','FontSize',16,'Location','Best')
end

% Particle Filter State Tracking
particle_track = [1, round(size(parameterStore,2)/2), size(parameterStore,2)];
colors = linspecer(length(particle_track));
% 1: Ksat, 2: Sy
subplot(4,1,4)
for jj = 1:length(particle_track)
    values = squeeze(parameterStore(:,particle_track(jj),1));
    plot(time_plot_days,parameterStore(:,particle_track(jj)),'Color',colors(jj,:),'LineWidth',2);
    set(gca,'FontSize',12)
    xlabel('Elapsed Time [days]','Interpreter','latex','FontSize',16);
    ylabel('$K_{sat}$~[$\mathrm{m \cdot s^{-1}}$]','Interpreter','latex','FontSize',16);
    hold on
end
legend('1-th','N/2-th','N-th','interpreter','latex','location','best');
exportgraphics(gcf,'Discharge_Ensemble.png','Resolution',300)
end

%% Initial Ensemble
% subplot(2,1,1)
% initial_depth = GW_table(col_int,row_in,1)
% plot(initial_depth,'linewidth',2,'color','black')
% xlabel('Particle','interpreter','latex','fontsize',12)
% ylabel('$h(0)$ [m/s]','interpreter','latex','fontsize',12)
% xlim([0,20])
% xlim([1,20])
% grid on
% hold on
% hold on
% plot([1:1:20],zeros(1,20),'linewidth',1.5,'linestyle','--')
% legend('Initial Parameter Ensemble','Inverse Problem Value','interpreter','latex')
% 
% 
% subplot(2,1,2)
% plot(parameterStore(1,:),'linewidth',2,'color','black')
% xlabel('Particle','interpreter','latex','fontsize',12)
% ylabel('$K_{sat}$ [m/s]','interpreter','latex','fontsize',12)
% xlim([0,20])
% xlim([1,20])
% grid on
% hold on
% hold on
% plot([1:1:20],ones(1,20)*(1/3600),'linewidth',1.5,'linestyle','--')
% legend('Initial Parameter Ensemble','Inverse Problem Value','interpreter','latex')

%% GW
% Recharge Values
TopBC_line = squeeze(TopBC(col_out,row_out,:));
line_plot(time_plot_days,'\mathrm{Elapsed~Time}','days',squeeze(GW(col_out,row_out,:)),'{h}_{\mathrm{GW}}','m',-TopBC_line*3600*1000,ChangeT/86400,'\bar{i}','\mathrm{mm\cdot h^{-1}}','GW Table at Chosen Cell')
exportgraphics(gcf,'Output/Figures/GW_Table_at_Outlet.png','Resolution','600')

%% Hydrograph Rate
line_plot(time_plot_days,'\mathrm{Elapsed~Time}','days',discharge*1000*3600*24,'Q','\mathrm{mm\cdot day^{-1}}',-TopBC_line*3600*1000*24,ChangeT/86400,'\bar{N}','\mathrm{mm\cdot day^{-1}}','Hydrograph')
exportgraphics(gcf,'Output/Figures/GW_Hydrograph_Rate.png','Resolution','600')
% line_plot(time_plot_days,1000*discharge*(nx*dx*ny*dy),'Q','\mathrm{L\cdot s^{-1}}')

%% Hydrograph
line_plot(time_plot_days,'\mathrm{Elapsed~Time}','days',discharge*1000*3600,'Q','\mathrm{mm\cdot h^{-1}}',-TopBC_line*3600*1000,ChangeT/86400,'\bar{N}','\mathrm{mm\cdot h^{-1}}','Hydrograph')
exportgraphics(gcf,'Output/Figures/GW_Hydrograph.png','Resolution','600')
% line_plot(time_plot_days,1000*discharge*(nx*dx*ny*dy),'Q','\mathrm{L\cdot s^{-1}}')
zzz = [GW_Storage, max(discharge'*1000*3600,0)];

%% Storage-Discharge
scatter_plot(1000*GW_Storage/(dx*dy*numel(X)),discharge*1000*3600*24,time_plot_days,'\bar{h}_{\mathrm{GW}}','Q','mm','mm \cdot day^{-1}','GW Storage-Discharge')
exportgraphics(gcf,'Output/Figures/GW_Storage_Discharge.png','Resolution','600')

%% Level Plots
nx_values = [1  round(nx/2)  nx];
ny_values = [1  round(ny/2)  ny];
% Hmin = min(min(min(1*z0)));
% Hmax = max(max(max(GW + 1*z0)));
% level_plots(time_plot_days,'Water Level Plots',GW,'z + P/\gamma',z0,'z',nx_values,ny_values,X,Y,Hmin,Hmax,'\mathrm{Head}','m')

%% GW Depth Plots
% Hmin = min(min(min(GW + 0*z0)));
% Hmax = max(max(max(GW + 0*z0)));
% level_plots(time_plot_days,'GW Depth Plots',GW,'P/\gamma',0*z0,'\mathrm{Distance~alongisde~hillslope}',nx_values,ny_values,X,Y,Hmin,Hmax,'\mathrm{Head Pressure}','m')

%% Normalized Plot
intervals = 3;
name_plot = 'Normalized_Storage';
figure(1)
obj = VideoWriter(['Output/Animations/' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
nx_values = [1  round(nx/2)  nx];
ny_values = [1  round(ny/2)  ny];
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end

for i = 1:(intervals+1)
    dx_time = length(time_plot_days)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    level_plot(time_plot_days(idx_time),100*GW(:,:,idx_time)./H,0*GW(:,:,idx_time),nx_values,ny_values,X,Y,'S/S_{{max}}','Ground',0,100,'S/S_{\mathrm{max}}','\%')
    f = getframe(gcf);
    writeVideo(obj,f);
end
obj.close();    
close all

%% GW Maps
intervals = 5;
rows = 2;
close all
for i = 1:(intervals+1)
    subplot(rows,ceil((intervals+1)/rows),i);
    dx_time = length(time_plot_days)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    surf_plot(max(max(max(GW))),time_plot_days(idx_time),'h_{GW}','m',GW(:,:,idx_time),1,1,256,0.9,0,[0,90],X,Y);
end
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
% %% 3D Field Groundwater Table
% field_plots(max(max(max(max(GW_table)))),time_plot_days,'Groundwater_Table_3D','GW','m',GW_table,Y(:),X(:),X,Y,1,1,0,[45,60])
exportgraphics(gcf,'Output/Figures/GW_Depths_Surfs.png','Resolution','600')
close all

%% Normalized Plots

nx_plot = [2  round(nx/2)  nx];
ny_plot = [2  round(ny/2)  ny];
intervals = 5;

profile_figures(intervals,time_plot_days,GW./H,X,Y,nx_plot,ny_plot,'Profile Figure','S/S_{\mathrm{max}}','-')


%% Normalized Plot
intervals = 15;
name_plot = 'Normalized_Storage';
figure(1)
obj = VideoWriter(['Output/Animations/' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
nx_values = [1  round(nx/2)  nx];
ny_values = [1  round(ny/2)  ny];
intervals = min(intervals,length(time_plot_days));
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end

for i = 1:(intervals+1)
    dx_time = length(time_plot_days)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    level_plot(time_plot_days(idx_time),100*GW(:,:,idx_time)./H,0*GW(:,:,idx_time)./H,nx_values,ny_values,X,Y,'S/S_{{max}}','Ground',0,100,'S/S_{\mathrm{max}}','\%')
    f = getframe(gcf);
    writeVideo(obj,f);
end
obj.close();    
close all

%% Field Groundwater
field_plots_3D(max(max(max(GW))),time_plot_days,'Groundwater_Field_2D','h_{GW}','m',GW,[0,90],'GW Depth [m]',X,Y)

%% 3D Field Groundwater Depth
field_plots_3D(max(max(max(GW))),time_plot_days,'Groundwater_Field_3D','h_{GW}','m',GW,[45,60],'GW Depth [m]',X,Y)
% %% 3D Fiedl Groundwater Table
% field_plots(max(max(max(max(GW_table)))),time_plot_days,'Groundwater_Table_3D','GW','m',GW_table,Y(:),X(:),X,Y,1,1,0,[45,60])

