% Particle Filter Framework
% Developer: Marcus Nobrega, Ph.D.
% Goal: For a given system dynamics, estimate states and outputs based on
% observations and modelling
% If your system is different, you need to change the parameters and states
%
% Parameters:
% Ksat,h: horizontal ksat [m/s]
% n: Van-Genutchen n [-]
% Poro: Effective porosity [-]
% Alpha: Van-Genutchen alpha [-] (negative quantity)
% irr_damping: A value between 0 and 1 such that irr = irr*irr_damping

% x = [Ksat_values(j), n_values(1), Poro_values(1), Alpha_values(1),, irr_damping(1)];

% Initializing the P.F
if time == 0

    % Read Particle Filter Inputs
    dir = obj.coupledModel.gwmodel.model_dir;
    input_data = readtable(dir,'Sheet','Assimilation_Data'); % Reading Input data from Excel

    PF_seed = table2array(input_data(1:6,2:5));

    % PF parameters
    parameter_std = table2array(input_data(1:6,3))';
    parameter_min = table2array(input_data(1:6,5))';
    parameter_max = table2array(input_data(1:6,6))';
    parameters = table2array(input_data(1:6,4))';
    initial_std_parameter = table2array(input_data(1:6,2))'; % Logical index showing which GW Depth gauges are used in the PF

    % Gauges Used in the PF
    gauges_used = logical(table2array(input_data(20:34,2)));

    dt = (ti - time); % time-step in seconds

    % Domain Discretization for the 2D Boussinesq Model
    nx = size(obj.coupledModel.gwmodel.parameter.Sy,2);
    ny = size(obj.coupledModel.gwmodel.parameter.Sy,1);

    % Parameter Definition
    % parameters = [obj.coupledModel.gwmodel.parameter.K(round(ny/2),round(nx/2)), ...
    %               obj.coupledModel.gwmodel.parameter.N(round(ny/2),round(nx/2)), ...
    %               obj.coupledModel.gwmodel.parameter.Poro(round(ny/2),round(nx/2)), ...
    %               obj.coupledModel.gwmodel.parameter.Alpha(round(ny/2),round(nx/2)), ...
    %               obj.coupledModel.gwmodel.parameter.Irr_d, ...
    %               obj.coupledModel.gwmodel.parameter.k(round(ny/2),round(nx/2))];

    % ---- Parameter Ranges and STD ---- %
    % parameter_std = [5*10^(-6), ...                           % GW ksat [m/s]
    %                  1*10^(-2), ...                           % n [-]
    %                  0*10^(-2), ...                           % porosity [-]
    %                  2*10^(-2), ...                           % alpha [1/m]
    %                  1*10^(-3), ...                         % irr_damping [-]
    %                  1*10^(-6)];                              % k [1/sec]
    % parameter_min = [1*10^(-7), ...                           % GW ksat [m/s]
    %                  1.00, ...                                % n [-]
    %                  0.385, ...                                % porosity [-]
    %                  -5.00, ...                               % alpha [1/m]
    %                  0.80, ...                                % irr_damping [-]
    %                  1*10^(-7)];                              % k [1/sec]
    % parameter_max = [1*10^(-2), ...                           % GW ksat [m/s]
    %                  5.00, ...                                % n [-]
    %                  0.385, ...                               % porosity [-]
    %                  -1, ...                                  % alpha [1/m]
    %                  1.2, ...                                 % irr_damping [-]
    %                  1*10^(-3)];                              % k [1/sec]

    % Initial Ensemble Std
    % initial_std_parameter = [5*10^(-5), ...                           % GW ksat [m/s]
    %                          10*10^(-2), ...                           % n [-]
    %                          0.0*10^(-2), ...                           % porosity [-]
    %                          20*10^(-2), ...                          % alpha [1/m]
    %                          2*10^(-1), ...                           % irr_damping [-]
    %                          0.9*10^(-5)];                              % k [1/sec]

    % initial_std_parameter(2) = 0;
    % initial_std_parameter(3) = 0;
    % initial_std_parameter(4) = 0;
    % initial_std_parameter(5) = 0;
    % initial_std_parameter(6) = 0;

    % parameter_std(2) = 0;
    % parameter_std(3) = 0;
    % parameter_std(4) = 0;
    % parameter_std(5) = 0;
    % parameter_std(6) = 0;

    %% System Description
    % Description of variables:
    % State = UZ Ksat, alpha, and n. GW Ksat
    % Output = Outflow discharge (mm / h) at the seepage face

    % Time-steps
    nt = 1; % Number of time-steps (this is one because we are only moving the model forward)

    % --- Process Std ---- %
    % GW heads are the first states
    x0_GW = obj.coupledModel.gwmodel.h(:)'; % [m]
    std_GW = [0.000]; % Process std of GW states

    % UZ moisture depth are the second states
    x0_S_UZ = obj.coupledModel.gwmodel.parameter.S_UZ(:)'; % [m]
    std_S_UZ = [0.000]; % Process std of GW states

    % Sy as the last state
    x0_Sy = obj.coupledModel.gwmodel.parameter.Sy(:)'; % [-]
    sdt_S_GW = [0.000]'; % Process std of Sy

    % Initial State
    x0 = [x0_GW, x0_S_UZ, x0_Sy];

    % --- Process Std
    % GW
    ne = length(x0_GW);
    idx_state(1,1) = ne;
    Std(1:ne,1) = std_GW(1);

    % UZ
    idx_state(1,2) = length(x0_GW) + length(x0_S_UZ);
    Std((ne+1):(idx_state(1,2)))= std_S_UZ(1);
    ne = idx_state(1,2) ;

    % Sy
    idx_state(1,3) = length(x0);
    Std((ne+1):(idx_state(1,3))) = sdt_S_GW(1);

    % Augmented State
    x0_aug = [x0,parameters];  % First entries are the concatenated states, last entries are the model parameters

    temp = ones(size(obj.coupledModel.gwmodel.z0));
    temp = temp(:);

    % Minumum and Maximum States
    state_min = [obj.coupledModel.gwmodel.z0(:);                                                 % bedrock elevation [m]
                0*temp;                                                                          % UZ storage [m]
                0*temp]';                                                                        % Sy [-]
    state_max = [obj.coupledModel.gwmodel.z0(:) + obj.coupledModel.gwmodel.parameter.H(:);       % bedrock + soil depth [m]
                10000*temp;                                                                      % No limit for storage [m]
                1*temp]';                                                                        % 1 for sy [m]

    %% Filter Parameters
    N = table2array(input_data(10,2));                    % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.
    flag_discharge = table2array(input_data(12,2));       % Indicates if we will use discharge as assimilated data
    n_outputs = flag_discharge + sum(gauges_used);
    process_std = Std;                                    % Process std for each state
    measurement_std_discharge = table2array(input_data(14,2));      % Measurement std. In this case, discharge in mm/h
    measurement_std_GW = table2array(input_data(16,2));   % Measurement std. In this case, GW depth in m
    measurement_std = [measurement_std_discharge,repmat(measurement_std_GW,1,n_outputs-1)];
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

    V = 0*(0.001^2); % define the variance of the initial state estimate

    % make the randomly generated particles from the initial prior gaussian distribution
    x_P = zeros(N,ns_aug);

    % % States required to save that are not in the augmented state
    % Sy_posterior = zeros(N,length(x0_GW)); % Specific yield

    % Indexes
    state_index = 1:length(x0); % From 1 to number of states
    parameter_index = ((length(x0)+1):1:ns_aug); % From number of states + 1 to states + parameters

    % State (Adding gaussian noise in the states, following V)
    % x_P(:,state_index) = x0_aug(:,state_index) + gaussian_noise(0,V,N,1); % Adding a gaussian noise in the input
    % Each single state with a different noise
    x_P(:,state_index) = x0_aug(:,state_index) + gaussian_noise(0,V,N,length(state_index)); % Adding a gaussian noise in the input

    % Parameters (Adding gaussian noise in the parameters, following
    % initial_std_parameter)
    % Each single state with a different noise
    % x_P(:,parameter_index) = x0_aug(:,parameter_index) + gaussian_noise(0,initial_std_parameter.^2,N,1); % Particle information related to the parameters
    x_P(:,parameter_index) = x0_aug(:,parameter_index) + gaussian_noise(0,initial_std_parameter.^2,N,length(parameter_index)); % Particle information related to the parameters

    %%%%% Manual Ensemble
    % factor = [0.8 0.9 1 1.1 1.2]';
    % x_P(:,parameter_index) = repmat(x0_aug(:,parameter_index),[N,1]).*repmat(factor,[1,length(parameter_index)]);

    % Applying a physical constraint to estimated states
    x_P(:,state_index) = max(x_P(:,state_index),state_min);
    x_P(:,state_index) = min(x_P(:,state_index),state_max); % Herein you can define what is a maximum value for the parameter

    % Applying a physical constraint in the parameters (augmented state)
    % x_P(:,parameter_index) = max(x_P(:,parameter_index),parameter_min); % Here you can enter physical meaning for the states or the bounds
    % x_P(:,parameter_index) = min(x_P(:,parameter_index),parameter_max); % Here you can enter physical meaning for the parameters or the bounds

    labels = {'Ksat [m/s]','n','Porosity','$\alpha$ [1/m]','$i_d$ [-]','$k$ [1/sec]'};
    % for ii = 1:length(initial_std_parameter)
    %     subplot(1,5,ii)
    %     plot(1:N,x_P(:,parameter_index(ii)),'linewidth',2,'color','black'); ylabel(labels{ii},'interpreter','latex'); xlabel('Particle index','interpreter','latex');
    % end
    % figure(4)
    plot_initial_ensemble(x_P(:,parameter_index),labels,10);
    exportgraphics(gcf,'Output/Figures/Initial_Ensemble.png','Resolution','600')
    savefig('Output/Figures/Initial_Ensemble.fig');
    % Mean Estimate of the State
    x0_est = mean(x_P(:,state_index));

    % Standard Deviation of the Estimated State
    std_begin_x = std(x_P(:,state_index));

    % Particle State
    x_P_update = x_P;

    %% Initial Values for the Simulation

    % Particle Output
    z_P_update = zeros(N,n_outputs); % Modeled output including discharge gw depth
    z_P_update_discharge = zeros(N,1);
    z_P_update_depth = zeros(N,n_outputs-1);

    % Particle Probability
    P_w = zeros(N,1); % Particle probability


end % Initialization is over

%% Loading Observed Output (z)
z_estimated = table2array(input_data(:,8:9)); % Time and Discharge
% Here you need to specify the observed data
obs_time = round(z_estimated(:,1)); % Time in minutes
discharge = z_estimated(:,2); % Discharge in mm/h
idx_time = find(ti/60 >= obs_time,1,'last'); % Finding the correct value for the current time (forward)
z_theoretical_discharge = discharge(idx_time); % mm/h

% Now we go for the GW depth
DX = obj.coupledModel.gwmodel.nx*obj.coupledModel.gwmodel.dx;
DY = obj.coupledModel.gwmodel.ny*obj.coupledModel.gwmodel.dy;
for j = 1:15 % Up to 15 points
    GW_Depth_Obs(:,j) = table2array(input_data(:,10 + (j-1)*3));
    x_well(j,1) = ceil(table2array(input_data(1,11 + (j-1)*3))/obj.coupledModel.gwmodel.dx); % Number of cells
    y_well(j,1) = ceil((DY - table2array(input_data(1,12 + (j-1)*3)))/obj.coupledModel.gwmodel.dy); % Number of cells
end
z_theoretical_depth = GW_Depth_Obs(idx_time,:);

% Taking away wells not used
z_theoretical_depth(~gauges_used) = [];
x_well(~gauges_used) = [];
y_well(~gauges_used) = [];
%% Loop through time
for t = 1:(nt)   % Since we have only one time, this only occurs once
    % z = z_theoretical_discharge + gaussian_noise(0,measurement_variance,1,1);
    z = [z_theoretical_discharge, z_theoretical_depth]; % No noise added into it
    %
    % % Physical Constraint in the observed value
    % z = max(z,0);
    % z = min(z,100000);

    % ------ Particle Filter ------ %
    for ii = 1:N % For all particles (This Could be Parallelizable)
        if time == 0
            % Original States
            obj.coupledModel.gwmodel.h = original_h_GW;
            obj.coupledModel.gwmodel.parameter.S_UZ = original_S_UZ;
            obj.coupledModel.gwmodel.parameter.Sy = original_Sy;
            % Generating an Ensemble for initial water depth
            % initial_depth(ii,1) = 1*rand();
            % obj.coupledModel.gwmodel.h(:) = ones(size(original_h_GW))*(initial_depth(ii,1)) + obj.coupledModel.gwmodel.z0;
        else
            % Initial states for the time-step
            current_states = states_posterior;
            obj.coupledModel.gwmodel.h = reshape(current_states(ii,(0+1):idx_state(1)),obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.nx);
            obj.coupledModel.gwmodel.parameter.S_UZ = reshape(current_states(ii,(idx_state(1)+1):idx_state(2)),obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.nx);
            obj.coupledModel.gwmodel.parameter.Sy = reshape(current_states(ii,(idx_state(2)+1):idx_state(3)),obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.nx);
        end

        % Original Irrigation
        idx = find(obj.coupledModel.gwmodel.boundaries.changeTQ<= time,1,'last');
        obj.coupledModel.gwmodel.boundaries.Q(:,:,idx) = original_irrigation(:,:,idx); % m/s

        % Particle Parameters
        parameters = x_P(ii,parameter_index); % New set of parameters

        % Assigning Parameters to the Model (GW)
        obj.coupledModel.gwmodel.parameter.K =  parameters(1)*ones(size(obj.coupledModel.gwmodel.parameter.K));
        obj.coupledModel.gwmodel.parameter.N =  parameters(2)*ones(size(obj.coupledModel.gwmodel.parameter.K));
        obj.coupledModel.gwmodel.parameter.Poro =  parameters(3)*ones(size(obj.coupledModel.gwmodel.parameter.K));
        obj.coupledModel.gwmodel.parameter.Alpha = parameters(4)*ones(size(obj.coupledModel.gwmodel.parameter.K));
        obj.coupledModel.gwmodel.parameter.Irr_d = parameters(5);

        % Make sure the boundaries have correct Ksat values
        obj.coupledModel.gwmodel.parameter.K(end,:) = 10000*parameters(1);


        % Run model forward
        big_number = 1e6;
        try
            obj.coupledModel.advance(time,ti);
            x_GW = obj.coupledModel.gwmodel.h(:)'; % [m]
            x_UZ = obj.coupledModel.gwmodel.parameter.S_UZ(:)'; % [m]
            x_Sy = obj.coupledModel.gwmodel.parameter.Sy(:)'; % [m]
            x_P_update(ii,state_index) = [x_GW, x_UZ, x_Sy]; % State
            x_P_update(ii,parameter_index) = parameters; % Parameters not changed
            % Estimated Output for each particle (discharge)
            if flag_discharge == 1
                z_P_update_discharge(ii,1) = max(1000*3600*obj.coupledModel.getDischarge,0); % mm/h
            else
                z_P_update_discharge(ii,1) = [];
            end
            % Estimated Output for each particle (GW depth)
            for kk = 1:size(z_theoretical_depth,2)
                z_P_update_depth(ii,kk) = obj.coupledModel.gwmodel.h(y_well(kk),x_well(kk)) - obj.coupledModel.gwmodel.z0(y_well(kk),x_well(kk)); % GW Depth in m
            end
            z_P_update(ii,:) = [z_P_update_discharge(ii,1), z_P_update_depth(ii,:)]; % Discharge and Depth
        catch ME
            x_P_update(ii,state_index) = nan*ones(1,ns); % States totally wrong
            x_P_update(ii,parameter_index) = parameters; % Parameters not changed
            % Estimated Output for each particle
            z_P_update(ii,:) = big_number*ones(1,n_outputs); % huge value for outflow (mm/h)
        end

        % -- Given the observation, what is the probability of the
        %    particle estimate?
        if t == 1
            % Initializing Weights
            P_w_0 = ones(N,1)*1/N; % all particles with the same weight
        else
            P_w_0 = P_w;
        end
        [P_w(ii,:)] = particle_probability(z_P_update(ii,:)',z',measurement_variance,flag_discharge,P_w_0(ii));
    end

    % Normalize to form a probability distribution (i.e. sum to 1).
    P_w = P_w./nansum(P_w);

    % Checking Particle Degeneracy
    if max(P_w) == 1
        % We have particle degeneracy, in this case we assume all particles now
        % have the same probability
        warning('Particle Degeneracy occuring')
        P_w = 1/N*ones(length(P_w),1);
    end


    parameter_prev = x_P(:,parameter_index);

    if isnan(nansum(P_w))
        warning('No particles satisfied the model constraints')
    end

    %% Resampling
    % From this new distribution, now we randomly sample from it to generate our new estimate particles
    % Now, let's consider a spin-up to start to update the particles
    time_spinup = 14; % Days;
    if time > time_spinup*1440*60 % Seconds
        time_span_calc = 0.5; % Hours
        dt_calc = time_span_calc/(dt/3600); 
        if mod(i, dt_calc) == 0
            [x_P,index] = resample_particles(P_w,x_P_update(:,state_index),process_variance,N,x_P_update(:,parameter_index),parameter_variance,z_P_update,z',measurement_std);
            P_w = ones(length(P_w),1)*1/N; % Reinitialize Weight
        end
    else
        x_P = x_P_update;
    end
    % clf
    % plot(x_P(:,parameter_index),'LineStyle','-','Color','red');
    % hold on
    % plot(parameter_prev,'LineStyle','-','Color','blue');
    % plot(1:1:N,1e-4*ones(1,N),'LineWidth',2,'LineStyle','--')
    % legend('Posterior','Prior')

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
    discharge_PF = z_P_update(:,1);
    GW_Depth_PF = z_P_update(:,2:end);
    Sy_PF = x_P(:,(idx_state(2)+1):idx_state(3));
    Sy_PF = nanmean(Sy_PF,2); % Saving only the average (> 2GB variables are not saved in matlab)
    parameter_PF = temp2;

    % Saving states for each particle for posterior
    states_posterior = x_P(:,state_index);
    parameters_posterior = x_P(:,parameter_index);

    % Plotting Particle Distribution
    if t == 1
        % close all

        % set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
    else
        clf
    end

    % tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    % nexttile
    % % Cell in which states will be plotted
    % cell_id = size(obj.coupledModel.gwmodel.z0)/2; % mid of the domain
    % state_id = 0 + cell_id(1) + (cell_id(2) - 1)*size(obj.coupledModel.gwmodel.z0,1);
    % particle_distribution(x_P(:,state_id) - obj.coupledModel.gwmodel.z0(cell_id(1),cell_id(2)),...
    %                         P_w,x_est(state_id) - obj.coupledModel.gwmodel.z0(cell_id(1),cell_id(2)),ti/60,'\mathrm{Probability}','\mathrm{Head~Pressure~at~half~of~the~domain}','-','m','Particles','Average State',[0,1])
    % hold off
    % % Kh, Sy
    % nexttile
    % % scatter_plot_scaled(x_P(:,parameter_index(2)),x_P(:,parameter_index(1)),'S_{\mathrm{y}}','K_{\mathrm{sat,h}}','-','m/s',P_w,strcat('t =~',num2str(round(ti/60,2)),'~min'))
    % % set(gca,'YScale','log')
    % % hold off
    %
    % % Plotting Discharge
    % nexttile
    % z_plot = z_P_update; z_plot(z_plot == big_number) = nan;
    % lim = [0, 1.5*max(max(z_P_update,z))];
    % if lim(2) == lim(1)
    %     lim(2) = lim(1) + 0.01;
    % end
    % particle_distribution(z_plot,...
    %                         P_w,z,ti/60,'\mathrm{Probability}','\mathrm{Discharge}','-','mm/h','Particle Output','Observed Output',lim)
    % pause(0.00001)

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

