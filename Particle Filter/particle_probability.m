function [P_w] = particle_probability(x,mean,variance,flag_discharge,P_w_0)
% Calculates the probability of a particle, given an observation
%
% Input:
% x: particle state
% mean: observed value
% variance: variance of the observations
% size: number
% Gaussian probability function
% Given an estimated output, what is the probability of the
% observation?

x_Z = variance;

n_wells = size(variance,1) - flag_discharge;

if flag_discharge == 1 && n_wells == 0
    alpha = 1;
elseif flag_discharge == 0 && n_wells > 0
   alpha = 0;
else
    alpha = 0.5; % We give alpha*100 percentage of importance for discharge
end

% Weights for each output
beta = 1 - alpha; % We give this importance for the average error of GW depth

% Gaussian Distribution
% P_w = (1./(sqrt(x_Z)*sqrt(2*pi))).* exp(-(x - mean).^2./(2.*x_Z));

% Method presented in "Estimation of Parameters in Groundwater Modeling by
% Particle Filter linked to the meshless local Petrov-Galerkin Numerical Method"
P_w = (P_w_0.* exp(-(x - mean).^2./(2.*x_Z)));

if flag_discharge == 0
    P_w_Q = 0;
else
    P_w_Q = alpha*P_w(1,1);
end

if n_wells > 0
    P_w_GW = 1/n_wells*beta.*P_w(2:end,1);
else
    P_w_GW = [];
end

% New Probability Vector
P_w = sum([P_w_Q; P_w_GW]);
end