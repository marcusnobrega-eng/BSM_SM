% Particle Filter Framework
% Developer: Marcus Nobrega, Ph.D.
% Goal: For a given system dynamics, estimate states and outputs based on
% observations and modelling


% Example: Catchment - Reservoir
% State: h - water depth [m]
% Output: Q - discharge [m3/s]
% Input: Inflow hydrograph [m3/s]

%% Tank Data
% Outflow = k1 * (h - h01)^k2 + k3 * (h - h02) ^ k4 
k1 = [1.35]';
k2 = [0.5]';
h01 = [0]';
k3 = [8]';
k4 = [1.5]';
h02 = [4]';
S0 = [0]'; % m3
A = [1000]'; % Tank area [m2]

parameters = [k1,k2,h01,k3,k4,h02]; % Model Parameters
parameter_std = [0.01 0.00 0.00 0.01 0.00 0.00]; % Model Parameter Std. 

% ----- Inflow Scenarios ----- %
dt = 15; % min
% Qin = [1 2 5 10 20 40 15 10 5 1]; % Inflow [m3/s]
for j = 1:length(A)
    % Qin = 20*ones(1,200);
    % Alternated Blocks + SCS Hydrograph
    [time_rainfall,i,~,~] = alternated_blocks(120,dt,1519,0.236,16,0.935,j*10,0);
    [~,Qin(:,j),~,~] = SCS_Hydrologic_Model(dt*60,720,i,dt,time_rainfall,[],[],[],2000,4000,0.2,90,0.02,0.035,0,0.015,0,1,1,0,50,0,0);
end
%% System Description
% Description of variables:
% State = Water Depth [m]
% Output = Outflow discharge (Qout - m3/s)
% Qout = k1 (h - h01)^k2 + k3(h - h02)^k3

% Time-steps
nt = length(Qin); % Number of time-steps
% Time Vector
time = 0:dt:nt; % min
% Initial Depth based on the initial volume
x0 = S0./A;

% Augmented State
x0_aug = [x0,parameters]; 
parameters_correct = parameters;
%% Filter Parameters
x = x0; % Initial actual state [m]
process_std = 0.01; % 
measurement_std = 0.01; % 
N = 100; % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.
no = [1]'; % Number of outputs
ns = [1]'; % Number of states

process_variance = [process_std.^2]'; % Noise variance in the system
measurement_variance = [measurement_std.^2]'; % Noise variance in the measurement
parameter_variance = [parameter_std.^2]; % Parameter variance used to estimate the gaussian noise
ns_aug = length(x0_aug); % Number of augmented states
%% Initializing Particles
% initilize our initial, prior particle distribution as a gaussian around
% the true initial value

V = 10*process_variance; % define the variance of the initial estimate

% make the randomly generated particles from the initial prior gaussian distribution

x_P = zeros(N,ns_aug);

% Indexes
state_index = 1:length(x0); % From 1 to number of states
parameter_index = ((length(x0)+1):1:ns_aug); % From number of states + 1 to states + parameters

% State
x_P(:,state_index) = x0_aug(state_index) + gaussian_noise(0,V,N,1); % Adding a gaussian noise in the input

% Parameters
x_P(:,parameter_index) = x0_aug(parameter_index) + gaussian_noise(0,100*parameter_variance,N,1); % Particle information related to the parameters

% Apply a physical constraint to estimated states
x_P(:,state_index) = max(x_P(:,state_index),0);
x_P(:,state_index) = min(x_P(:,state_index),1000); % Herein you can define what is a maximum value for the parameter

% Apply a physical constraint in the parameters (augmented state)
x_P(:,parameter_index) = max(x_P(:,parameter_index),1e-3); % Here you can enter physical meaning for the states or the bounds
x_P(:,parameter_index) = min(x_P(:,parameter_index),1000); % Here you can enter physical meaning for the parameters or the bounds

% Mean Estimate of the State
x0_est = mean(x_P(:,state_index));

% Standard Deviation of the Estimated State
std_begin_x = std(x_P(:,state_index));

% Particle State
x_P_update = x_P;

% Standard deviation of output
z0 = output_dynamics(parameters,x_P(:,state_index)); % Output calculated with estimated states from the particles
std_begin_z = std(z0); % Estimated particle standard deviation on ouput
z0_est = mean(z0); % Mean Estimated particle output

%% Initial Values for the Simulation
% ---- Initial Values --- %
Inflow = Qin(1); % Inflow hydrograph from a correct (could be uncertain) hydrological model [m3/s]

% Preallocating arrays
x_expected = zeros(ns,nt);        % State
z_expected = zeros(ns,nt);        % Output
x_estimated  = zeros(ns,nt);      % Estimatated deepth
z_estimated = zeros(ns,nt);       % Estimated discharge
std_estimated_x = zeros(ns,nt);   % Std deviation of states
std_estimated_z = zeros(ns,nt);   % Std deviation of discharge
x_theoretical = zeros(ns,nt);     % Theoretical
z_theoretical = zeros(ns,nt);     % Theoretical

% State
xe0 = x0; % Initial Estimate of x

% Output
z0 = output_dynamics(parameters,x0); % Output estimated with x0

% Looping parameters
x_previous = x0; 
z_previous = z0;
Inflow_previous = Inflow;

% Saved values
x_expected(1,1) = x0;
z_expected(1,1) = z0;
x_estimated(1,1) = x0_est;
z_estimated(1,1) = z0_est;
std_estimated_x(1,1) = std_begin_x;
std_estimated_z(1,1) = std_begin_z;

x_perfect(1,1) = x0;
z_perfect(1,1) = z0;
x_previous_perfect = x_previous;
z_previous_perfect = z0;
%% Loop through time
for t = 1:(nt)
    Inflow = Qin(t);

    % No Uncertainy Case (Perfect State and Perfect Output)
    x_theoretical = dynamic_model(x_previous_perfect,[Inflow_previous;z_previous_perfect], ...
                [A],dt);
    z_theoretical = output_dynamics(parameters_correct,x_theoretical); % Calculated with perfect state and with correct parameters

    % State
    % x = max(tank_dynamics(Inflow_previous,z_previous,A,dt,x_previous) + ...
        % 0*gaussian_noise(0,process_variance,1),0); % State + Noise + Physical Constraint

    x = dynamic_model(x_previous,[Inflow_previous;z_previous], ...
                [A],dt);

    % Output: theoretical value + gaussian noise 
    z = z_theoretical + gaussian_noise(0,measurement_variance,1,1);

    % Physical Constraint in the observed value
    z = max(z,0);
    z = min(z,100000);

    % z = max(outflow_discharge(k1,k2,h01,k3,k4,h02,x) + ...
        % 0*gaussian_noise(0,measurement_variance,1),0); % Output + Noise + Physical Constraint

    % ------ Particle Filter ------ %
    x_P_update(:,state_index) = dynamic_model(x_P(:,state_index),[Inflow_previous;z_previous], ...
                                              [A],dt);
    % x_P_update = max(tank_dynamics(Inflow_previous,z_previous,A,dt,x_P),0); % State + Noise + Physical Constraint

    % Estimated Output for each particle
    parameters = x_P(:,parameter_index); % New set of parameters   
    % z_P_update = outflow_discharge(k1,k2,h01,k3,k4,h02,x_P_update); % No noise
    z_P_update = output_dynamics(x_P(:,parameter_index),x_P_update(:,state_index)); % No noise

    % -- Given the observation, what is the probability of the
    %    particle estimate?     
    [P_w] = particle_probability(z_P_update',z,measurement_variance);
    
    % Normalize to form a probability distribution (i.e. sum to 1).
    P_w = P_w./sum(P_w);
    
    %% Resampling 
    % From this new distribution, now we randomly sample from it to generate our new estimate particles

    [x_P] = resample_particles(P_w,x_P_update(:,state_index),process_variance,N,parameters,parameter_variance);   

    % Apply a physical constraint to estimated states
    x_P(:,state_index) = max(x_P(:,state_index),0);
    x_P(:,state_index) = min(x_P(:,state_index),10000);
    
    % Apply a physical constraint in the parameters (augmented state)
    x_P(:,parameter_index) = max(x_P(:,parameter_index),1e-3); % Here you can enter physical meaning for the states or the bounds
    x_P(:,parameter_index) = min(x_P(:,parameter_index),1000); % Here you can enter physical meaning for the parameters or the bounds
    
    % The final estimate is some metric of the final resampling, such as
    % the mean value or std deviation

    x_est = mean(x_P(:,state_index));  % Estimated mean state 
    z_est = mean(z_P_update); % Estimated mean output
    std_est_x = std(x_P(:,state_index)); std_est_z = std(z_P_update); % Estimated ouput std
    
    % Plotting Particle Distribution
    if t == 1
        close all

        set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
    else
        clf
    end
    subplot(1,2,1)
    particle_distribution(x_P_update(:,state_index),P_w,x,dt*(t-1),'\mathrm{Probability}','X_p','-','m')
    subplot(1,2,2)    
    parameter_distribution(parameters,parameters_correct,dt*(t-1),'\mathrm{Parameter}','\mathrm{Value}','-','-')
    pause(0.0001);

    % Save data in arrays for later plotting
    x_expected(1,t) = x;
    z_expected(1,t) = z;
    x_estimated(1,t) = x_est;
    z_estimated(1,t) = z_est;
    std_estimated_x(1,t) = std_est_x;
    std_estimated_z(1,t) = std_est_z;

    % Theoretical case
    x_perfect(1,t) = x_theoretical;
    z_perfect(1,t) = z_theoretical;

    % Refreshing Inputs
    x_previous = x; 
    z_previous = z;
    Inflow_previous = Inflow;
    x_previous_perfect = x_theoretical;
    z_previous_perfect = z_theoretical;
    
end

%% Calculate Fitness Metrics
[r2, ia, nse, kge, PBIAS, rmse, mae] = fitness_metrics(x_expected,x_estimated);
metrics = [r2, ia, nse, kge, PBIAS, rmse, mae];
% Define the headers
headers = {'r2', 'ia', 'nse', 'kge', 'PBIAS', 'rmse', 'mae'};

% Create a table with the metrics
Table = table(metrics(1), metrics(2), metrics(3), metrics(4), metrics(5), metrics(6), metrics(7), ...
          'VariableNames', headers);

% Display the table in MATLAB
disp(Table);
disp('Metrics calculated based on observed and mean modeled value.')

% Export the table to a CSV file
writetable(Table, 'performance_metrics.csv');

%% Plotting
close all
[Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps();
t = [0:1:(nt-1)]*dt;
figure(1);
clf
plot(t, x_perfect, '-','linewidth',2,'Color',pallete.green_colors(1,:));
hold on
plot(t, x_expected, '-','linewidth',2,'Color',pallete.red_colors(1,:));
hold on
% errorbar(t, x_estimated, std_estimated_x, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
plot(t, x_estimated, '-','linewidth',2,'Color',pallete.blue_colors(2,:));
hold on
set(gca,'FontSize',12); set(gcf,'Color','White');
fill([t, fliplr(t)], [x_estimated + std_estimated_x, fliplr(x_estimated - std_estimated_x)], pallete.blue_colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
xlabel('Elapsed time [min]','Interpreter','latex','FontSize',16); ylabel('Depth [m]','Interpreter','latex','FontSize',16)
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Perfect Model','Observed', 'Particle filter estimate','Interpreter','latex','FontSize',16,'Location','Best')
exportgraphics(gcf,'Output.png','Resolution',300)

figure(2)
clf
plot(t, Qin, '-','linewidth',3,'Color',pallete.blue_colors(1,:));
% plot(t,x_estimated,'-.','linewidth',2,'Color',pallete.red_colors(2,:));
hold on
plot(t, z_perfect, '-.','linewidth',3,'Color',pallete.green_colors(1,:));
hold on
plot(t, z_expected, '-.','linewidth',3,'Color',pallete.red_colors(1,:));
hold on
plot(t, z_estimated, ':','linewidth',3,'Color',pallete.blue_colors(2,:));

% Define the upper and lower bounds for shading
fill([t, fliplr(t)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(2,:), ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Elapsed time [min]','Interpreter','latex','FontSize',16); ylabel('Discharge [$\mathrm{m^3 \cdot s^{-1}}$]','Interpreter','latex','FontSize',16)
% Customize ticks
set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
legend('Inflow','Perfect $Q_{\mathrm{out}}$','Observed $Q_{\mathrm{out}}$', 'Particle Filter $Q_{\mathrm{out}}$','Interpreter','latex','FontSize',16,'Location','Best')
exportgraphics(gcf,'State.png','Resolution',300)

