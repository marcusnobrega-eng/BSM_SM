% Optimization Function for the Optimization Problem
% Semi-distributed Model Coupled with Channel Routing
% Developer: Marcus Nobrega, Ph.D.
% 4/24/2024

function [Obj_fun,modeled_discharge,GW_Depth_Modeled_Total] = Opt_function(x, dir, flag_opt_fun,flag_GW)
% Add Paths
addpath Setups/
addpath GW/
addpath Output/
addpath Extra/

cmf = General_Setup_Calibrator(x,dir); % Calling 
try
%% Observed Wells 
% x_well = [12 12 12 12];
% y_well = [50 41 28 20];

%%
[discharge,GWdepth,dx,dy,DX,DY,nx,ny] = run_boussinesqModel(cmf);

%% Observed GW
input_data = readtable(dir,'Sheet','Assimilation_Data'); % Reading Input data from Excel
% Now we go for the GW depth
for j = 1:15 % Up to 15 points
    GW_Depth_Obs(:,j) = table2array(input_data(:,10 + (j-1)*3));
    x_well(j,1) = ceil(table2array(input_data(1,11 + (j-1)*3))/dx); % Number of cells
    y_well(j,1) = ceil((DY - table2array(input_data(1,12 + (j-1)*3)))/dy); % Number of cells
end

GW_Depth_Obs(1,:) = []; % First value not used in the calibration

% Gauges Used in the PF
gauges_used = logical(table2array(input_data(20:34,2)));

% Taking away wells not used
GW_Depth_Obs(:,~gauges_used) = [];
x_well(~gauges_used) = [];
y_well(~gauges_used) = [];


%% Modeled GW
for i = 1:length(x_well)
    data = reshape(GWdepth,ny,nx,size(GWdepth,2));
    GW_Depth_Modeled(:,i) = squeeze(data(y_well(i),x_well(i),:));
end

%% Read Observed Outlet Hydrograph
input_table = readtable('Observed_Hydrograph.csv');
Q_obs_table = table2array(input_table(:,:));
time_Q_obs = Q_obs_table(:,1);
Q_obs = Q_obs_table(:,2:end); % only last entries

%% Extract only modeled flows at observed time-steps
t_step_begin = 673; % Number of steps to start
dt = cmf.coupledModel.gwmodel.dt;
time_begin = t_step_begin*dt/60; % min
t_step_begin_obs = find(time_Q_obs <= time_begin,1,'last'); 

% Cutting Observed Values
Q_obs(1:t_step_begin_obs,:) = [];
time_Q_obs(1:t_step_begin_obs,:) = [];

% Cutting Modeled Values
modeled_discharge = discharge;
discharge = discharge';
discharge(1:t_step_begin,:) = [];

GW_Depth_Modeled_Total = GW_Depth_Modeled;

GW_Depth_Modeled(1:t_step_begin,:) = [];
GW_Depth_Obs(1:t_step_begin,:) = [];

% % Extracting Only Values at Observations
% steps_Qobs = ceil(time_Q_obs*60/dt - t_step_begin + 1); % steps where we have observed data
% discharge = discharge(steps_Qobs,1);

Q_m = (1)*1000*3600*discharge;

flag_GW = 1;
%% Objective functions
if flag_GW ~= 1
    [Obj_fun] = metrics(Q_obs,Q_m,flag_opt_fun);
else
    [Obj_function_Flow] = metrics(Q_obs,Q_m,flag_opt_fun);
    for i = 1:size(GW_Depth_Obs,2)
        [Obj_function_GW(i)] = metrics(GW_Depth_Obs(:,i),GW_Depth_Modeled(:,i),flag_opt_fun);
    end

    Obj_fun = 1/2*(Obj_function_Flow + mean(Obj_function_GW));
end


catch ME
    Obj_fun = inf;
    modeled_discharge = nan;
    GW_Depth_Modeled_Total = nan;
end
end

function [Obj_fun] = metrics(Q_obs,Q_m,flag_opt_fun)
if flag_opt_fun == 1 % NSE
    Qobs_avg = mean(Q_obs);
    % Nash-Suctclife-Efficiency and Objective function
    NSE = 1 - sum((Q_obs - Q_m).^2)/(sum((Q_obs - Qobs_avg).^2)); % We want to maximize it
    % Vol Error
    % Vol_error = (sum(Q_m*dt) - sum(Q_obs*dt))/(sum(Q_obs*dt));
    % Obj_fun = -0.8*NSE -0.2*Vol_error; % Therefore, we want to minimize it
    % Backing to NSE
    Obj_fun = -NSE;
end

if flag_opt_fun == 2 % RMSE
    n_elements = length(Q_obs);
    RMSE = sum((Q_obs - Q_m).^2/n_elements);
    Obj_fun = RMSE; % minimize RMSE
end

if flag_opt_fun == 3 % R2
    r2 = corrcoef(Q_obs,Q_m);
    Obj_fun = -r2(1,2); % maximize r2
end

if flag_opt_fun == 4 % Peak Flow
    z1 = max(Q_obs);
    z2 = max(Q_m);
    Obj_fun = abs(z1 - z2); % Absolute value of peak flows
end

end













