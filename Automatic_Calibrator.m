%%%% Optimization Problem - Main Script
% Solving Optimization Problem to determine hydrologic soil parameters
% Developer: Marcus Nóbrega, Ph.D.
% Date: 04/25/2024
% Goal: Use optimization techniques to calibrate a coupled 1DUZ 2DGW model

%% Theory behind the method

%%% Length of Decision Vector
% Given the UZ parameters of:
% ksat: saturated hydraulic conductivity
% porosity: total porosity 
% alpha: Van-Genuthchen parameter
% n: Van-Genutchen parameters
% and assuming that ksat(UZ) = ksat(GW), 
% we define the decision vector
% x = [ksat, porosity, alpha, n]^T that we want to find
% such that we minimize the error between the modeled and 
% observed hydrograph
% x = [x1, x2 ... xnp], where np is the number of parameters

%%% Minimum and Maximum values for decision variables (e.g., constraints)
% Given a identity matrix I of order np (I_np) and given a collumn
% minimum (xmin) and maximum (xmax) vectors, we can write the constraints
% as a matrix multiplication, such that:
% I_np*xmin <= x <= I_np*xmax

%%% Optimization function
% Let's assume that the stocastic or deterministic optimization problem
% gives us a ramdom or specific decision vector x. The optimization
% function depends of values values of x, such that we can say that:
% Opt_function = f(x)
% We probably won't be able to write an explicit equation for f(x), so that
% we will write an algorithm to obtain Opt_function.
% This algorithm will work as follows: (a) first choose an decision vector
% x, (b) then run the coupled model and retrieve the outlet
% hydrograph, (c) then choose parameters to match outlet flows with the
% observed flows, and finally (d) do that by minimizing 
% an error function.

%%% Optimization Problem
% minimize Opt_function = - NSE
% subject to:
% I_ndv * xmin <= x <= I_ndv * xmax

%% Input data
% Defining min and max constraints
xmin = [1e-6, ...                % Ksat [m/s]
        1, ...                   % N [-]
        0.37, ...                % Porosity [-]
        -5, ...                  % Alpha [1/m]
        1e-6, ...                % k [1/sec]
        0.90];                   % Irrigation damping [-]
xmax = [5e-3 ...                 % Ksat [m/s]
        5, ...                   % N [-]
        0.37, ...                % Porosity [-]
        -1, ...                  % Alpha [1/m]
        5e-4, ...                % k [1/sec]
        1];                      % Irrigation damping [-]

%% Checking errors
warning('off')
if length(xmin) ~= length(xmax)
    error('Dimension of xmin and xmax are different. Please correct it.')
end
%% Number of Decision Variables
n_dv = length(xmax); % Number of decision variables
I_dv = ones(n_dv,n_dv); % Identity matrix

%% Genetic Algorithm Optimization Problem
pop_size = 100;
generations = 20;
stall_generations = 10;
tolerance = 1e-6;
stall_limit = 24*60*60; % If no improvement is found in this time in seconds, the algorithm stops.
nvars = n_dv;
A = []; % Inequality constraint given by A*x <= b
b = []; 
Aeq = []; % Equality constraint given by Aeq*x = beq
beq = []; % 
lb = xmin;
ub = xmax;
nonlcon = []; % Non-linear constraints of x (function of x)
intcon = []; % Integer only constraints of x
% options = optimoptions('ga','PlotFcns', {@gaplotbestf,@gaplotstopping},'PopulationSize',pop_size,'Generations',generations,'StallTimeLimit',stall_limit,'MaxStallGenerations',stall_generations,'FunctionTolerance',tolerance);
options = optimoptions('ga','PopulationSize',pop_size,'Generations',generations,'StallTimeLimit',stall_limit,'MaxStallGenerations',stall_generations,'FunctionTolerance',tolerance,'UseParallel',true);
flag_opt_fun = 1;
% Excel Input
dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Input_Data_GW_Model_PERTH.xlsx';
fun = @(x)Opt_function(x, dir,flag_opt_fun);
% GA problem
[x,fval,exitflag,output,population,scores] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);

%% Run model with optimized parameters and Plot Results
% % Run Model
% Optimized_Model
% 
% % Retrieve Outlet Hydrograph and Observed Hydrograph
% %%%%% Outlet hydrograph
% flow_outlet = Q(:,outlet_pos); % Ony flow at the outlet;
% 
% %%%%%  Read Observed Outlet Hydrograph
% input_table = readtable('Observed_Hydrograph.csv');
% Q_obs_table = table2array(input_table(:,:));
% time_Q_obs = Q_obs_table(:,1);
% Q_obs = Q_obs_table(:,2:end); % only last entries
% time_step_Q_obs = (time_Q_obs(2,1) - time_Q_obs(1,1)); % minutes
% 
% %%%%%  Extract only modeled flows at observed time-steps
% steps_Qobs = ceil(time_Q_obs*60/dt + 1); % steps where we have observed data
% Q_m = flow_outlet(steps_Qobs);
% % Plot Results
% plot_optimal_results


