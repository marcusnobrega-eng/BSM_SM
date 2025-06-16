% Recharge Delay through a UZ reservoir
% Developer: Marcus Nobrega, Ph.D.
% 
% Input:
% I: irrigation [m/s]
% N: recharge rate
% k: damping factor [1/sec]
% dt: time-step [sec]
% S0: UZ storage [mm] at time t
% ST: UZ storage [mm] at time t + dt
%
% Dynamics:
% ds/dt = I - N
% ST - S0 = dt(I - N)
% ST = S0 + dt(I - N)
% N = K * S0
function [N,ST] = irrigation_delay(I,k,dt,S0)
    % ds/dt = I - N
    % ST - S0 = dt(I - N)
    % ST = S0 + dt(I - N)
    % N = K * S0
    if max(max(I)) ~= 0
        ttt = 1;
    end
    N = -k.*S0; % m/s
    ST = S0 + dt*(-I + N);

end