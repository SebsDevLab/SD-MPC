%% File Name: LPV_2_Tank.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: Non-linear cascaded two tank system modeled as LPV system. 
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
% 
% Outputs:
%   sys: The controllable system with its describing matrices and further
%   features
%
% Notes: 
% 


function sys = LPV_2_Tank()
    sys.name = "LPV_2_Tank";
    %% Nonlinear two tank state-space Model from [1]
    A_theta = @(theta_1, theta_2) [-0.904 * theta_1, 0; 
                                    0.904 * theta_1, -0.508 * theta_2];
    B = [0.258;
         0];
    C = [0,1];
    D = [0];
    
    sys.model.A = A_theta;
    sys.model.B = B;
    sys.model.C = C;
    sys.model.D = D;
    
    %% Calculate Sampling-Time for a given linearization point 
    x_1 = 1;    % Water level tank 1
    x_2 = 1;    % Water level tank 2
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);
    A = A_theta(theta_1, theta_2);

    sys_c = ss(A, B, C, D);
    poles = eig(sys_c);     % Calculate poles of the continous system
    tau = 1./poles;
    tau_fast = min(abs(tau));   % fastest pole
    
    f_tau = 1/(2*pi*tau_fast); 
    f_samp = 10*f_tau;

    T_samp = 1/f_samp;
    sys.T_samp = T_samp;

    %% Assume linearized State-space Model for parametric model 
    x_1 = 5;
    x_2 = 15;
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);

    A = A_theta(theta_1, theta_2);
    sys2 = ss(A, B, C, D);
    sys_d = c2d(sys2, T_samp);

    sys.param_model.A_M = sys_d.A;
    sys.param_model.B_M = sys_d.B;
    sys.param_model.C_M = sys_d.C;
    sys.param_model.D_M = sys_d.D;

    %% LPV state-space model of complete system (for simulation) 
    % sim_model_1: general senario: system data-collection == system simulation
    % sim_model_2: robust senario: system data-collection != system simulation

    sys.sim_model_1.A_sim = sys.model.A;
    sys.sim_model_1.B_sim = sys.model.B;
    sys.sim_model_1.C_sim = sys.model.C;
    sys.sim_model_1.D_sim = sys.model.D;
    
    sys.sim_model_2.A_sim = @(theta_1,theta_2) 0.9 * A_theta(theta_1, theta_2);
    sys.sim_model_2.B_sim = sys.model.B;
    sys.sim_model_2.C_sim = sys.model.C;
    sys.sim_model_2.D_sim = sys.model.D;

    %% System constants
    sys.nx = size(sys.model.A(1,1), 1); % Number of states of the system
    sys.nu = size(sys.model.B, 2);      % Number of inputs of the system
    sys.ny = size(sys.model.C, 1);      % Number of outputs of the system

    sys.nx_M = size(sys.param_model.A_M, 1); % Number of states of the parametric model
    sys.nu_M = size(sys.param_model.B_M, 2); % Number of inputs of the parametric model
    sys.ny_M = size(sys.param_model.C_M, 1); % Number of outputs of the parametric model

    
    
    %% Constraints
    sys.constraints.u_min = 0;
    sys.constraints.u_max = 22;
    sys.constraints.x_min = 0;
    sys.constraints.x_max = 100;
    sys.constraints.y_min = 0;
    sys.constraints.y_max = 100;
    
    
    %% Initial Condition
   
    
    
    
end