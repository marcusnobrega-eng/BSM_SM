function [x_P,index] = resample_particles(P_w,x_P_update,x_S,N,parameters,parameter_variance,z_PF,z,measurement_std)
% Input:
%
% P_w: particle probability
% x_P_update: particle state
% x_S: process variance
% N: number of particles
% parameters: parameters used in the model
% parameter_variance: variance in parameters
% measurement_std: measurement std of all outputs


% Avoid Ensemble Deterioration Phenomenon "An iterative particle filter
% approach for coupled hydro-geophysical inversion of a controlled infiltration experiment"

index = zeros(1,N); % Index showing which particles will be resampled

Neff = 1/(sum(P_w.^2));
[Neff, length(P_w)/2]
if Neff < 2*(length(P_w)/2)
    parameters_update = 0*parameters;

    rn = 1/N*rand; % Random number between 0 and 1/N
    cumprob = cumsum(P_w); % Cumulative probability of Pw
    x_P = 0*x_P_update; % Initializing particle states

    flag_resample_method = 1;
    % 1: Systematic Resampling
    % 2: Multinomial Resampling
    % 3: Residual resampling
    try
        for i = 1:N
            if flag_resample_method == 1
                % Systematic Resampling
                sampled_number = rn + (i-1)*1/N;
                index(1,i) = find(sampled_number <= cumprob,1,'first');
                x_P(i,:) = x_P_update(index(1,i),:); % No noise in States
                parameters_update(i,:) = parameters(index(1,i),:); % Parameters Replicated
            elseif flag_resample_method == 2
                % Draw a random number
                sampled_number = rand;
                index(1,i) = find(sampled_number <= cumprob,1,'first');
                x_P(i,:) = x_P_update(index(1,i),:); % No noise in States
                parameters_update(i,:) = parameters(index(1,i),:); % Parameters Replicated
            elseif flag_resample_method == 3
                % Draw a random number
                sampled_number = rn + (i-1)*1/N;
                if P_w(i) > 1/N
                    index(1,i) = i;
                else
                    index(1,i) = find(sampled_number <= cumprob,1,'first');
                end
                x_P(i,:) = x_P_update(index(1,i),:); % No noise in States
                parameters_update(i,:) = parameters(index(1,i),:); % Parameters Replicated
            end
        end
    catch
        ttt = 1;
    end
    % Find Repeated Parameters
    [~, ~, rowIdx] = unique(parameters_update, 'rows');
    counts = accumarray(rowIdx, 1);
    if N <= 3
        error('Please use more than 3 particles')
    end
    threshold_particle_min = max(ceil(N/10),2);
    threshold_particle_min = 2; %%%%%%
    threshold_particle_max = N;
    logical_vector = (counts(rowIdx) > threshold_particle_min) & (counts(rowIdx) < threshold_particle_max); % Index of particles that had parameters repeated more than once

    for ii = 1:length(counts)
        temp = rowIdx == ii;
        count_id = sum(temp);
        begin = find(rowIdx == ii,1,'first');
        for jj = 1:count_id
            logical_vector(begin + (jj-1),1) = 0;
            if jj > threshold_particle_min && jj < threshold_particle_max
                logical_vector(begin + (jj-1),1) = 1;
            end
        end
    end

    % Add noise only in cells that were resampled more than once
    x_P = x_P + logical_vector.*gaussian_noise(0,x_S,size(x_P,1),size(x_P,2));

    if ~isempty(parameters_update)

        % Add noise only in parameters that had states resampled more than once
        %%%% Method 1
        % parameters_update_std = std(parameters_update);
        % Lambda = 0.5; % Inflation Factor
        % parameter_variance_resample = parameters_update_std*Lambda;

        %%% Method 2
        parameter_update_std = 1.5*sqrt(var(parameters,P_w));
        parameters_update = parameters_update + 1*logical_vector(:,1).*gaussian_noise(0,parameter_update_std.^2,size(parameters,1),size(parameters,2));
        % parameters = parameters + 1*logical_vector(:,1).*gaussian_noise(0,parameter_variance,size(parameters,1),size(parameters,2));

        % % Check wether modeled outputs mean are close enough to \sqrt(std +
        % % error)
        % z_PF_mean = nanmean(z_PF)';
        % z_PF_std = nanstd(z_PF)';
        % PF_capacity = sqrt(z_PF_std.^2 + (measurement_std.^2)');
        % absolute_error = abs(z_PF_mean - z);
        % PF_check = PF_capacity./absolute_error;
        % % Add noise only in parameters that had states resampled more than once
        % parameters_update = parameters_update + 1*logical_vector(:,1).*gaussian_noise(0,parameter_variance,size(parameters,1),size(parameters,2));
        % % parameters = parameters + 1*logical_vector(:,1).*gaussian_noise(0,parameter_variance,size(parameters,1),size(parameters,2));
        % % Particle Inflation Method
        % if min(min(PF_check)) < 1
        %     idx_error = find(PF_check == min(min(PF_check)),1,'first');
        %     z_particles = z_PF(:,idx_error);
        %     idx_particles = double(abs(z_particles - z(idx_error)) > absolute_error(idx_error)); % Particles with error larger than tolerated error
        %     parameter_variance = 0*parameter_variance; % Adding 10-fold variance
        %     parameters_update = parameters_update + idx_particles.*gaussian_noise(0,parameter_variance,size(parameters,1),size(parameters,2));
        % end
        x_P = [x_P parameters_update];
    else
        parameters_update = [];
        x_P = x_P;
    end

else
    x_P = [x_P_update parameters];
end
end