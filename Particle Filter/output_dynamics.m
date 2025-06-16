function [Output] = output_dynamics(parameters,previous_state) % Output
% Calculates the model output
%
% Input:
% parameters = [parameter1, parameter2, ... parameterN] collects all output
% parameters
% previous_state: previous model state
%
% Example: Tank outflow devices
%
% Equations:
% Qout = Cd Aef sqrt(2g(h - h01)) + Cds Lef (h - h02) ^ (3/2)
% Qout = k1 (h - h01)^ k2 + k2 (h - h02)^k3
% parameters [k1;k2,h01;k3;k4;h02]
Output = parameters(:,1).*(max(previous_state - parameters(:,3),0)).^(parameters(:,2)) + ...
         parameters(:,3).*(max(previous_state - parameters(:,6),0)).^parameters(:,5);
end