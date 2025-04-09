%% File Name: Sinus.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: A sinus as reference trajectory swining around the value y_0
% Sources: 
%
%
% Inputs:
% T_sim: Simulation time
% ny: number of outputs 
% T_fut: Prediction horizon for a sufficient number of steps over the simulation horizon
% ini_len: Number of steps on y_0 to fill u_past and y_past of the data-driven component
% sin_period_len: length of the sinus period
% y_0: initial value of reference trajectory
% y_sin: amplitude of the sinus
%
%
% Outputs:
%   ref: the reference trajectory
%
% Notes: 
% 
function ref = Sinus(T_sim, ny, T_fut, ini_len, ~, ~, sin_period_len, y_0, ~, y_sin)
    ref = y_0 * ones(ny, T_sim + T_fut);
    lin = linspace(0, 2*pi, sin_period_len);
    ref2 = [zeros(1, ini_len), lin];
    for i=1:round(T_sim/sin_period_len)
        ref2 = [ref2, lin(2:end)];
    end
    ref2 = y_sin*sin(ref2);
    ref = ref + ref2(1:length(ref));
end