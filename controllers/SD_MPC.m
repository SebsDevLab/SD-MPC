%% File Name: SD_MPC.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: Controller class for SD-MPC
% Sources: 
% [1] - Sebastian Zieglmeier, et.al., "Semi-Data-Driven Model Pparamictive
%       Control: A Physics-Informed Data-Driven Control Approach", 
%       https://doi.org/10.48550/arXiv.2504.00746
%
% Properties: 
% T_d: Number of data collection points
% T_ini: Number of initialization points
% T_fut: Prediction horizon
% R, Q: cost matrices for input and output
% lambda_ini, lambda_g: regularization hyperparameter
% H_u, H_y: Hankel matrices for input and output data
% sys: system with parametric model 
% 
% Inputs Step function: 
% obj: Object
% x0: inital state for parametric model
% u_past, y_past: past T_ini data values for input and output
% y_ref: reference trajectory for prediction horizon
%
% Outputs Step function:
% u_fut: computed control input sequence
% x_M_value: future predicted states of the parametric model
% y_M_value: future predicted outputs of the parametric model
% yalmip_error_flag: error flag for optimization errors
% 
%  
% Notes: 
% 


classdef SD_MPC < handle
    properties
        % General control variables:
        T_d
        T_ini
        T_fut
        T
        R
        Q
        lambda_ini
        lambda_g

        % Data driven component
        H_u
        H_y
        
        % Parametric Model 
        A
        B
        C
        D
        
        % System dimensions
        nx
        nu
        ny
        nx_M
        nu_M
        ny_M

        % System boundaries
        u_min
        u_max
        x_min
        x_max
        y_min
        y_max
    end
    
    methods
        function obj = SD_MPC(T_d, T_ini, T_fut, R, Q, lambda_ini, lambda_g, H_u, H_y, sys)
            % General control variables
            obj.T_d = T_d;
            obj.T_ini = T_ini;
            obj.T_fut = T_fut;
            obj.T = T_ini + T_fut;
            obj.R = R;
            obj.Q = Q;
            obj.lambda_ini = lambda_ini;
            obj.lambda_g = lambda_g;
            
            % Data driven component
            obj.H_u = H_u;
            obj.H_y = H_y;
            
            % Parametric model 
            obj.A = sys.param_model.A_M;
            obj.B = sys.param_model.B_M;
            obj.C = sys.param_model.C_M;
            obj.D = sys.param_model.D_M;

            obj.nx = sys.nx; % Number of states
            obj.nu = sys.nu; % Number of inputs
            obj.ny = sys.ny; % Number of outputs
            
            obj.nx_M = sys.nx_M; % Number of states of parametric model
            obj.nu_M = sys.nu_M; % Number of inputs of parametric model
            obj.ny_M = sys.ny_M; % Number of outputs of parametric model

            obj.u_min = sys.constraints.u_min * ones(T_fut, sys.nu);
            obj.u_max = sys.constraints.u_max * ones(T_fut, sys.nu);
            obj.x_min = sys.constraints.x_min * ones(T_fut+1, sys.nx);
            obj.x_max = sys.constraints.x_max * ones(T_fut+1, sys.nx);
            obj.y_min = sys.constraints.y_min * ones(T_fut, sys.ny);
            obj.y_max = sys.constraints.y_max * ones(T_fut, sys.ny);
            
        end
        


        function [u_fut, x_M_value, y_M_value, yalmip_error_flag] = step(obj, x0, u_past, y_past, y_ref)
            verbose = false; 

            %% Initialize Optimization Variables
            
            % Data driven component
            g = sdpvar(obj.T_d - obj.T + 1, 1);
            y_D = sdpvar(obj.T_fut, obj.ny);
            u_D = sdpvar(obj.T_fut, obj.nu);
            sigma = sdpvar(obj.T_ini, obj.ny);

            % Parametric model 
            y_M = sdpvar(obj.T_fut, obj.ny_M);
            x_M = sdpvar(obj.T_fut+1, obj.nx_M);
            u_M = sdpvar(obj.T_fut, obj.nu);
            
            % Combined variables  
            u = sdpvar(obj.T_fut, obj.nu);
            y = sdpvar(obj.T_fut, obj.ny);

            % Reshape Hankel matrizes H_u and H_y
            U_p = obj.H_u(1:obj.nu * obj.T_ini, :);
            U_f = obj.H_u(obj.nu * obj.T_ini + 1:end, :);
            Y_p = obj.H_y(1:obj.ny * obj.T_ini, :);
            Y_f = obj.H_y(obj.ny * obj.T_ini + 1:end, :);

            %% Cost function
            cost = obj.Q*(y - y_ref)' * (y - y_ref) + obj.R * (u' * u) + obj.lambda_ini * (sigma' * sigma) + obj.lambda_g * norm(g, 1);
            % cost = obj.Q*(y(2:end,1) - y_ref(2:end,1))' * (y(2:end,1) - y_ref(2:end,1)) + obj.R * (u' * u) + obj.lambda_ini * (sigma' * sigma) + obj.lambda_g * norm(g, 1);
            % cost = obj.Q*(y - y_ref)' * (y - y_ref) + obj.R * (u' * u) + obj.lambda_ini * norm(sigma,1) + obj.lambda_g * norm(g, 1);

            %% Constraints
            % Parametric model: 
            constraints = [x_M(1,:)' == x0];
            for k =1:1:obj.T_fut
                constraints = [constraints, x_M(k+1,:)' == obj.A*x_M(k,:)' + obj.B*u_M(k,:)'];
                constraints = [constraints, y_M(k,:)' == obj.C*x_M(k,:)' + obj.D*u_M(k,:)'];
            end

            % Data driven component: 
            constraints = [constraints;
                U_p * g == u_past;
                U_f * g == u_D;
                Y_p * g == y_past + sigma;
                Y_f * g == y_D 
            ];

            % General control and constraints:
            constraints = [constraints;
                y == y_D + y_M;
                u == u_D;
                u == u_M;
                y(end,1) == y_ref(end,1); % Terminal constraint
                u >= obj.u_min;
                u <= obj.u_max;
                y >= obj.y_min;
                y <= obj.y_max;
            ];


            %% Solve the problem
            options = sdpsettings('verbose', verbose, 'solver', 'mosek');
            diagnostics = optimize(constraints, cost, options);
            yalmip_error_flag = 0;
            if diagnostics.problem ~= 0
                % warning('YALMIP had an issue with the optimization. Problem Code: %d', diagnostics.problem);
                yalmip_error_flag = 1;
            end
            u_fut = value(u);
            y_M_value = value(y_M);
            x_M_value = value(x_M);
        end
    end
end
