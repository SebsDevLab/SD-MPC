%% File Name: discretize_LPV.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: Function for discretizing the LPV system for simulating it on a given operating point. 
% Sources: 
% [1] - Teppa-Garrán, et.al,, "Liquid level tracking for a coupled tank system 
%       using quasi–lpv control.", Ingenius 33 (2025): 15-26.
% [2] - Sebastian Zieglmeier, et.al., "Semi-Data-Driven Model Pparamictive
%       Control: A Physics-Informed Data-Driven Control Approach", 
%       https://doi.org/10.48550/arXiv.2504.00746
%
%
%
% Inputs:
% A_theta, B_c, C_c, D_c: system properties of continous LPV system
% x: state of the system for linearizing the LPV system at the given state
% for simulating it
% T_samp: Sampling time for discretization
% Outputs:
%   A, B, C, D: discrete system dynamic matrices
%
% Notes: 
% 


function [A, B, C, D] = discretize_LPV(A_theta, B_c, C_c, D_c, x, T_samp)
    % Ensuring no division by zero
    x_1 = max(x(1,1), 1);
    x_2 = max(x(1,2), 1);
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);
    % Build discrete system
    A = A_theta(theta_1, theta_2);
    sys2 = ss(A, B_c, C_c, D_c);
    sys_d = c2d(sys2, T_samp);
    A = sys_d.A;
    B = sys_d.B; 
    C = sys_d.C; 
    D = sys_d.D;
end
    