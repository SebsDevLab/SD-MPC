%% File Name: Smooth_Step.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: A smooth step as reference trajectory 
% Sources: 
%
%
% Inputs:
% T_sim: Simulation time
% ny: number of outputs 
% T_fut: Prediction horizon for a sufficient number of steps over the simulation horizon
% ini_len: Number of steps on y_0 to fill u_past and y_past of the data-driven component
% smooth_len: Smoothening of the step e.g. to satisfy system limitations
% y_0: initial value of reference trajectory
% y_step: final value of the step
%
% Outputs:
%   ref: the reference trajectory
%
% Notes: 
% 
function ref = Smooth_Step(T_sim, ny, T_fut, ini_len, smooth_len, ~, ~, y_0, y_step, ~)
    ref = y_0 * ones(ny, T_sim + T_fut);
    step_end = y_step-(y_step-y_0)*0.1;
    ref(ny, ini_len-1:ini_len+smooth_len-3) = linspace(y_0, step_end, smooth_len+2-3);
    ref(ny, ini_len+smooth_len-2) = y_step-(y_step-y_0)*0.05;
    ref(ny, ini_len+smooth_len-1) = y_step-(y_step-y_0)*0.02;
    ref(ny, ini_len+smooth_len) = y_step-(y_step-y_0)*0.01;
    ref(ny, ini_len+smooth_len+1:T_sim + T_fut) = ones(size(ref(1, ini_len+smooth_len+1:T_sim + T_fut)))*y_step;
end