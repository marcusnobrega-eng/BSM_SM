% Particle Filter Framework
% Developer: Marcus Nobrega, Ph.D.
% Goal: For a given system dynamics, estimate states and outputs based on
% observations and modelling
% If your system is different, you need to change the parameters and states


% Initializing the P.F
if time == 0
dt = (ti - time); % time-step in seconds

% Parameter Definition
parameters = [obj.coupledModel.uzmodel.parameter.K(1,1), ...
              obj.coupledModel.uzmodel.parameter.Alpha(1,1), ...
              obj.coupledModel.uzmodel.parameter.N(1,1), ...
              obj.coupledModel.gwmodel.parameter.K(1,1)];

% ---- Parameter Ranges and STD ---- %
parameter_std = [1*10^(-6), ...                           % UZ ksat [m/s] 
                 0.2, ...                                 % UZ Alpha [1/m]
                 0.2, ...                                 % UZ n [-]
                 1*10^(-6), ...                           % GW ksat [m/s]
                ];

parameter_min = [1e-6, ...                                  % UZ ksat [m/s] 
                 2.33, ...                                  % UZ Alpha [1/m]
                 1.5, ...                                   % UZ n [-]
                 1e-6, ...                                  % GW ksat [m/s]
                ];

parameter_max = [1e-4, ...                                  % UZ ksat [m/s] 
                 4.33,   ...                                % UZ Alpha [1/m]
                 5.0, ...                                   % UZ n [-]
                 1e-4, ...                                  % GW ksat [m/s]
                ];

%% System Description
% Description of variables:
% State = UZ Ksat, alpha, and n. GW Ksat
% Output = Outflow discharge (mm / h) at the seepage face

% Time-steps
nt = 1; % Number of time-steps (this is one because we are only moving the model forward)

% Soil moisture content in UZ / matric head
x0_UZ = [obj.coupledModel.uzmodel(:).WC'];
std_UZ = [0.000]; % Process std of UZ states

% --- Process Std ---- %
nb = 1; % First state (begin)
ne = length(obj.coupledModel.uzmodel(:).WC'); % Number of states for the first state (end)
idx_state(1,1) = ne; % Index representing the end of the first states
Std(nb:ne,1) = std_UZ(1); % Assigning standard deviation to the first parameter

% GW heads are the states
x0_GW = obj.coupledModel.gwmodel.h(:)'; % [m] (second state)
std_GW = [0.000]; % Process std of GW states

% --- Process Std
nb = ne; 
ne = ne + length(obj.coupledModel.gwmodel.h(:)');
idx_state(2,1) = ne;
Std(nb:ne,1) = std_GW(1);

% Initial State
x0 = [x0_UZ, x0_GW];

% Augmented State
x0_aug = [x0,parameters];  % First entries are the concatenated states, last entries are the model parameters

nz = length(obj.coupledModel.uzmodel(:).h');

% Minumum and Maximum States
state_min = [(0.0001)*ones(1,nz)'; ...           % Minimum Water Content [-]
             obj.coupledModel.gwmodel.z0(:)]';   % bedrock elevation [m]

state_max = [(1)*ones(1,nz)'; ...                % Maximum Moisture Content [-]
             obj.coupledModel.gwmodel.z0(:) + obj.coupledModel.uzmodel.H(:)]'; % bedrock + soil depth [m]

%% Filter Parameters
N = 300;                     % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.
process_std = Std;           % Process std for each state
measurement_std = 0.05;      % Measurement std. In this case, discharge in mm/h

x = x0; % Initial actual state [m]
process_variance = [process_std.^2]';               % Noise variance in the system
measurement_variance = [measurement_std.^2]';       % Noise variance in the measurement
parameter_variance = [parameter_std.^2];            % Parameter variance used to estimate the gaussian noise
ns = length(x0);                                    % Number of states
no = 1;                                             % Number of outputs 
ns_aug = length(x0_aug);                            % Number of augmented states

%% Initializing Particles
% Initilize our initial, prior particle distribution as a gaussian around
% the true initial value

V = 0*(0.001^2); % define the variance of the initial estimate

initial_std_parameter = [1e-05, ...                                 % UZ ksat [m/s] 
                         0.2, ...                                   % UZ Alpha [1/m]
                         0.2, ...                                   % UZ n [-]
                         1e-5, ...                                  % GW ksat [m/s]
                        ];

% make the randomly generated particles from the initial prior gaussian distribution
x_P = zeros(N,ns_aug);

% States required to save that are not in the augmented state
Sy_posterior = zeros(N,length(x0_GW)); % Specific yield 
h_UZ_posterior = zeros(N,length(x0_UZ));

% Indexes
state_index = 1:length(x0); % From 1 to number of states
parameter_index = ((length(x0)+1):1:ns_aug); % From number of states + 1 to states + parameters

% State (Adding gaussian noise in the states, following V)
x_P(:,state_index) = x0_aug(state_index) + gaussian_noise(0,V,N,1); % Adding a gaussian noise in the input

% Parameters (Adding gaussian noise in the parameters, following
% initial_std_parameter)
x_P(:,parameter_index) = x0_aug(parameter_index) + gaussian_noise(0,initial_std_parameter.^2,N,1); % Particle information related to the parameters

% Applying a physical constraint to estimated states
x_P(:,state_index) = max(x_P(:,state_index),state_min);
x_P(:,state_index) = min(x_P(:,state_index),state_max); % Herein you can define what is a maximum value for the parameter

% Applying a physical constraint in the parameters (augmented state)
x_P(:,parameter_index) = max(x_P(:,parameter_index),parameter_min); % Here you can enter physical meaning for the states or the bounds
x_P(:,parameter_index) = min(x_P(:,parameter_index),parameter_max); % Here you can enter physical meaning for the parameters or the bounds

% Mean Estimate of the State
x0_est = mean(x_P(:,state_index));

% Standard Deviation of the Estimated State
std_begin_x = std(x_P(:,state_index));

% Particle State
x_P_update = x_P;

%% Initial Values for the Simulation

% Particle Output
z_P_update = zeros(N,1); % Modeled output

% Particle Probability
P_w = zeros(N,1); % Particle probability

% Domain Discretization for the 2D Boussinesq Model
nx = size(obj.coupledModel.gwmodel.parameter.Sy,2);
ny = size(obj.coupledModel.gwmodel.parameter.Sy,1);
end % Initialization is over

%% Loading Observed Output (z)
load z_estimated.csv
obs_time = round(z_estimated(:,1)); % Time in minutes
discharge = z_estimated(:,2); % Discharge in mm/h
idx_time = find(time/60 >= obs_time,1,'last'); % Finding the correct value for the current time
z_theoretical = discharge(idx_time); % mm/h

%% Loop through time
for t = 1:(nt)   % Since we have only one time, this only occurs once 
    % z = z_theoretical + gaussian_noise(0,measurement_variance,1,1);
    z = z_theoretical; % No noise added into it

    % Physical Constraint in the observed value
    z = max(z,0);
    z = min(z,100000);

    % ------ Particle Filter ------ %
    for ii = 1:N % For all particles (This Could be Parallelizable)        
        if time == 0
            % Original States
            obj.coupledModel.uzmodel(:).WC = original_WC_UZ';
            obj.coupledModel.gwmodel.h(:) = original_h_GW;
            obj.coupledModel.gwmodel.parameter.Sy = original_Sy;
            obj.coupledModel.uzmodel(:).h = original_h_UZ';
        else
            % Initial states for the time-step
            current_states = states_posterior; 
            obj.coupledModel.uzmodel(:).WC = current_states(ii,1:idx_state(1))';
            obj.coupledModel.gwmodel.h = reshape(current_states(ii,(idx_state(1)+1):idx_state(2)),obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.nx);
            obj.coupledModel.gwmodel.parameter.Sy = reshape(Sy_posterior(ii,:),obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.nx);
            obj.coupledModel.uzmodel(:).h = h_UZ_posterior(ii,:)';
        end
         
        % Particle Parameters
        parameters = x_P(ii,parameter_index); % New set of parameters   

        % Assigning Parameters to the Model (UZ)
        obj.coupledModel.uzmodel.parameter.K = parameters(1)*ones(size(obj.coupledModel.uzmodel.parameter.K));
        obj.coupledModel.uzmodel.parameter.Alpha = parameters(2)*ones(size(obj.coupledModel.uzmodel.parameter.Alpha));
        obj.coupledModel.uzmodel.parameter.N = parameters(3)*ones(size(obj.coupledModel.uzmodel.parameter.N));

        % Assigning Parameters to the Model (GW)
        obj.coupledModel.gwmodel.parameter.K = parameters(4)*ones(size(obj.coupledModel.gwmodel.parameter.K));

        % Make sure the boundaries have correct Ksat values
        obj.coupledModel.gwmodel.parameter.K(end,2:(end-1)) = 1000*parameters(4);

        % Run model forward
        big_number = 1e6;
        try
            obj.coupledModel.advance(time,ti);
            x_UZ = [obj.coupledModel.uzmodel(:).WC'];
            x_GW = obj.coupledModel.gwmodel.h(:)'; % [m]
            x_P_update(ii,state_index) = [x_UZ, x_GW]; % State            
            % Estimated Output for each particle
            z_P_update(ii,1) = 1000*3600*obj.coupledModel.getDischarge; % mm/h 
            % Saving state-dependent Sy
            Sy_posterior(ii,:) = obj.coupledModel.gwmodel.parameter.Sy(:)'; % 
            % Saving h UZ
            h_UZ_posterior(ii,:) = obj.coupledModel.uzmodel.h(:)'; % 
        catch ME
            x_P_update(ii,state_index) = nan*ones(1,ns); % States totally wrong
            % Estimated Output for each particle
            z_P_update(ii,1) = big_number; % huge value for outflow (mm/h)
            Sy_posterior(ii,:) = obj.coupledModel.gwmodel.parameter.Sy(:)'; % 
            % Saving h UZ
            h_UZ_posterior(ii,:) = obj.coupledModel.uzmodel.h(:)'; % 
        end

        % -- Given the observation, what is the probability of the
        %    particle estimate?     
        [P_w(ii,1)] = particle_probability(z_P_update(ii,1)',z,measurement_variance);        
    end   

    
    % Normalize to form a probability distribution (i.e. sum to 1).
    P_w = P_w./nansum(P_w);

    if isnan(nansum(P_w))
        warning('No particles satisfied the model constraints')
    end
    
    %% Resampling 
    % From this new distribution, now we randomly sample from it to generate our new estimate particles

    [x_P,index] = resample_particles(P_w,x_P_update(:,state_index),process_variance,N,x_P_update(:,parameter_index),parameter_variance);   

    % Consider the Sy and h from index
    Sy = Sy_posterior;
    h_UZ = h_UZ_posterior;
    for jj = 1:N
        Sy_posterior(jj,:) = Sy(index(jj),:); % Using resampled info 
        h_UZ_posterior(jj,:) = h_UZ(index(jj),:); % Using resampled info 
    end

    % Apply a physical constraint to estimated states
    x_P(:,state_index) = max(x_P(:,state_index),state_min);
    x_P(:,state_index) = min(x_P(:,state_index),state_max);
    
    % Apply a physical constraint in the parameters (augmented state)
    x_P(:,parameter_index) = max(x_P(:,parameter_index),parameter_min); % Here you can enter physical meaning for the states or the bounds
    x_P(:,parameter_index) = min(x_P(:,parameter_index),parameter_max); % Here you can enter physical meaning for the parameters or the bounds
        
    % The final estimate is some metric of the final resampling, such as
    % the mean value or std deviation
    temp = x_P(:,state_index); temp(isnan(nan)) = nan;
    temp1 = z_P_update; temp1(temp1 == big_number) = nan;
    temp2 = x_P(:,parameter_index); temp2(temp2 == 0) = nan;

    x_est = nanmean(temp);  % Estimated mean state 
    z_est = nanmean(temp1); % Estimated mean output
    parameter_est = nanmean(temp2); % Average Parameter Estimation
    std_est_x = nanstd(temp1); 
    std_est_z = nanstd(temp2); % Estimated ouput std

    % Particle filter outputs for saving
    discharge_PF = z_P_update;
    parameter_PF = temp2;

    % Saving states for each particle for posterior
    states_posterior = x_P(:,state_index);
    parameters_posterior = x_P(:,parameter_index);
    Sy_posterior = Sy_posterior(index,:); % Using resampled info 
    
    % Plotting Particle Distribution
    if t == 1
        % close all

        % set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
    else
        clf
    end    

    tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');                    
    nexttile
    % Cell in which states will be plotted
    cell_id = size(obj.coupledModel.gwmodel.z0)/2; % mid of the domain
    state_id = idx_state(1) + cell_id(1) + (cell_id(2) - 1)*size(obj.coupledModel.gwmodel.z0,1);
    particle_distribution(x_P(:,state_id) - obj.coupledModel.gwmodel.z0(cell_id(1),cell_id(2)),...
                            P_w,x_est(state_id) - obj.coupledModel.gwmodel.z0(cell_id(1),cell_id(2)),ti/60,'\mathrm{Probability}','\mathrm{Head~Pressure~at~half~of~the~domain}','-','m','Particles','Average State',[0,1])
    
    % Plotting WC at the surface
    nexttile    
    state_id = idx_state(1);
    lim = [min(min(x_P(:,state_id))),1.0*max(max(x_P(:,state_id)))];
    if lim(2) == lim(1)
        lim(2) = lim(1) + 0.01;
    end
    particle_distribution(x_P(:,state_id),...
                            P_w,x_est(state_id),ti/60,'\mathrm{Probability}','\mathrm{WC~at~the~surface}','-','-','Particles','Average State',lim)     
    
    % Plotting Water Content
    nexttile
    state_id = round(idx_state(1)/2);
    lim = [min(min(x_P(:,state_id))),1.0*max(max(x_P(:,state_id)))];
    if lim(2) == lim(1)
        lim(2) = lim(1) + 0.01;
    end
    particle_distribution(x_P(:,state_id),...
                            P_w,x_est(state_id),ti/60,'\mathrm{Probability}','\mathrm{WC~at~half~of~the~UZ~col.}','-','m','Particles','Average State',lim) 

    % Plotting Discharge
    nexttile
    z_plot = z_P_update; z_plot(z_plot == big_number) = nan;
    lim = [0, 1.5*max(max(x_P(:,state_id)),z)];
    if lim(2) == lim(1)
        lim(2) = lim(1) + 0.01;
    end
    particle_distribution(z_plot,...
                            P_w,z,ti/60,'\mathrm{Probability}','\mathrm{Discharge}','-','mm/h','Particle Output','Observed Output',lim) 
    pause(0.00001)

    % Refreshing Inputs
    x_previous = x; 
    z_previous = z;    
end

% %% Calculate Fitness Metrics
% [r2, ia, nse, kge, PBIAS, rmse, mae] = fitness_metrics(x_expected,x_estimated);
% metrics = [r2, ia, nse, kge, PBIAS, rmse, mae];
% % Define the headers
% headers = {'r2', 'ia', 'nse', 'kge', 'PBIAS', 'rmse', 'mae'};
% 
% % Create a table with the metrics
% Table = table(metrics(1), metrics(2), metrics(3), metrics(4), metrics(5), metrics(6), metrics(7), ...
%           'VariableNames', headers);
% 
% % Display the table in MATLAB
% disp(Table);
% disp('Metrics calculated based on observed and mean modeled value.')
% 
% % Export the table to a CSV file
% writetable(Table, 'performance_metrics.csv');

%% Plotting
% close all
% [Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps();
% t = [0:1:(nt-1)]*dt;
% figure(1);
% clf
% plot(t, x_perfect, '-','linewidth',2,'Color',pallete.green_colors(1,:));
% hold on
% plot(t, x_expected, '-','linewidth',2,'Color',pallete.red_colors(1,:));
% hold on
% % errorbar(t, x_estimated, std_estimated_x, 'o-', 'LineWidth', 1, 'Color',pallete.red_colors(1,:), 'MarkerSize', 4, 'DisplayName', 'Particle Estimate');
% plot(t, x_estimated, '-','linewidth',2,'Color',pallete.blue_colors(2,:));
% hold on
% set(gca,'FontSize',12); set(gcf,'Color','White');
% fill([t, fliplr(t)], [x_estimated + std_estimated_x, fliplr(x_estimated - std_estimated_x)], pallete.blue_colors(3,:), ...
%     'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
% xlabel('Elapsed time [min]','Interpreter','latex','FontSize',16); ylabel('Depth [m]','Interpreter','latex','FontSize',16)
% % Customize ticks
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
% legend('Perfect Model','Observed', 'Particle filter estimate','Interpreter','latex','FontSize',16,'Location','Best')
% exportgraphics(gcf,'Output.png','Resolution',300)
% 
% figure(2)
% clf
% plot(t, Qin, '-','linewidth',3,'Color',pallete.blue_colors(1,:));
% % plot(t,x_estimated,'-.','linewidth',2,'Color',pallete.red_colors(2,:));
% hold on
% plot(t, z_perfect, '-.','linewidth',3,'Color',pallete.green_colors(1,:));
% hold on
% plot(t, z_expected, '-.','linewidth',3,'Color',pallete.red_colors(1,:));
% hold on
% plot(t, z_estimated, ':','linewidth',3,'Color',pallete.blue_colors(2,:));
% 
% % Define the upper and lower bounds for shading
% fill([t, fliplr(t)], [z_estimated + std_estimated_z, fliplr(z_estimated - std_estimated_z)], pallete.blue_colors(2,:), ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Error Range of Time Series 1');
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Elapsed time [min]','Interpreter','latex','FontSize',16); ylabel('Discharge [$\mathrm{m^3 \cdot s^{-1}}$]','Interpreter','latex','FontSize',16)
% % Customize ticks
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Garamond');
% legend('Inflow','Perfect $Q_{\mathrm{out}}$','Observed $Q_{\mathrm{out}}$', 'Particle Filter $Q_{\mathrm{out}}$','Interpreter','latex','FontSize',16,'Location','Best')
% exportgraphics(gcf,'State.png','Resolution',300)

