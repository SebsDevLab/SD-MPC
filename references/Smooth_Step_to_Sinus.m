%% File Name: Smooth_Step_to_Sinus.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: A smooth step with an overlying sinus as reference trajectory 
% Sources: 
%
%
% Inputs:
% T_sim: Simulation time
% ny: number of outputs 
% T_fut: Prediction horizon for a sufficient number of steps over the simulation horizon
% ini_len: Number of steps on y_0 to fill u_past and y_past of the data-driven component
% smooth_len: Smoothening of the step e.g. to satisfy system limitations
% step_len: Length of the Step without overlying sinus 
% sin_period_len: length of the sinus period
% y_0: initial value of reference trajectory
% y_step: final value of the step
% y_sin: amplitude of the sinus
%
%
% Outputs:
%   ref: the reference trajectory
%
% Notes: 
% 

function ref = Smooth_Step_to_Sinus(T_sim, ny, T_fut, ini_len, smooth_len, step_len, sin_period_len, y_0, y_step, y_sin)
    % Smooth Step:
    ref = y_0 * ones(ny, T_sim + T_fut);
    step_end = y_step-(y_step-y_0)*0.1;
    ref(ny, ini_len-1:ini_len+smooth_len-3) = linspace(y_0, step_end, smooth_len+2-3);
    ref(ny, ini_len+smooth_len-2) = y_step-(y_step-y_0)*0.05;
    ref(ny, ini_len+smooth_len-1) = y_step-(y_step-y_0)*0.02;
    ref(ny, ini_len+smooth_len) = y_step-(y_step-y_0)*0.01;
    ref(ny, ini_len+smooth_len+1:T_sim + T_fut) = ones(size(ref(1, ini_len+smooth_len+1:T_sim + T_fut)))*y_step;
    
    % Sinus:
    lin = linspace(0, 2*pi, sin_period_len);
    ref2 = [zeros(1, ini_len+step_len), lin];
    for i=1:round(T_sim/sin_period_len)
        ref2 = [ref2, lin(2:end)];
    end
    ref2 = y_sin*sin(ref2);
    ref = ref + ref2(1:length(ref));
end