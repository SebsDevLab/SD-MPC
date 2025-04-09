%% File Name: Step.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: A step with as reference trajectory 
% Sources: 
%
%
% Inputs:
% T_sim: Simulation time
% ny: number of outputs 
% T_fut: Prediction horizon for a sufficient number of steps over the simulation horizon
% ini_len: Number of steps on y_0 to fill u_past and y_past of the data-driven component
% y_0: initial value of reference trajectory
% y_step: final value of the step
%
%
% Outputs:
%   ref: the reference trajectory
%
% Notes: 
% 
function ref = Step(T_sim, ny, T_fut, ini_len, ~, ~, ~, y_0, y_step, ~)
    ref = y_0 * ones(ny, T_sim + T_fut);
    ref(ny, ini_len+1:T_sim + T_fut) = ones(size(ref(1, ini_len+1:T_sim + T_fut)))*y_step;
end