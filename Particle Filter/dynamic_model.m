function [state] = dynamic_model(previous_state,input,parameters,time_step)
% Calculates the state at k = k + 1
% Input data:
%
% x_p: previous state
% input: collects all model inputs (can be dynamic)
% parameters: collect all model parameters
% time_step: model time-step
%
% Exemple: Tank Model
% Parameters = [A] - Tank Area
% Input = [Inflow ; Outflow] - m3/s
% Dynamical Equation
% dh/dt = 1/A (Qin - Qout)
% h(t + dt) = h(t) + dt/A(Qin - Qout)
state = previous_state + time_step/parameters(1)*(input(1) - input(2));
end