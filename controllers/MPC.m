%% File Name: MPC.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: Controller class for SD-MPC
% Sources: 
% [1] - Sebastian Zieglmeier, et.al., "Semi-Data-Driven Model Pparamictive
%       Control: A Physics-Informed Data-Driven Control Approach", 
%       https://doi.org/10.48550/arXiv.2504.00746
%
% Properties: 
% T_fut: Prediction horizon
% R, Q: cost matrices for input and output
% sys: system with parametric model 
% 
% Inputs Step function: 
% obj: Object
% x0: inital state for parametric model
% y_ref: reference trajectory for prediction horizon
%
% Outputs Step function:
% u_fut: computed control input sequence
% x_M_value: dummy vector for using same step function framework as SD-MPC
% y_M_value: dummy vector for using same step function framework as SD-MPC
% yalmip_error_flag: error flag for optimization errors
% 
%  
% Notes: 
% 


classdef MPC < handle
    properties
        % General control variables:
        T_fut
        R
        Q
        
        % Parametric Model 
        A
        B
        C
        D
        
        % System dimensions
        nx
        nu
        ny

        % System boundaries
        u_min
        u_max
        x_min
        x_max
        y_min
        y_max
    end
    
    methods
        function obj = MPC(~, ~, T_fut, R, Q, ~, ~, ~, ~, sys)
            % General control variables
            obj.T_fut = T_fut;
            obj.R = R;
            obj.Q = Q;
            
            % Parametric model 
            obj.A = sys.param_model.A_M;
            obj.B = sys.param_model.B_M;
            obj.C = sys.param_model.C_M;
            obj.D = sys.param_model.D_M;

            obj.nx = sys.nx_M; % Number of states
            obj.nu = sys.nu_M; % Number of inputs
            obj.ny = sys.ny_M; % Number of outputs

            obj.u_min = sys.constraints.u_min * ones(T_fut, sys.nu);
            obj.u_max = sys.constraints.u_max * ones(T_fut, sys.nu);
            obj.x_min = sys.constraints.x_min * ones(T_fut+1, sys.nx);
            obj.x_max = sys.constraints.x_max * ones(T_fut+1, sys.nx);
            obj.y_min = sys.constraints.y_min * ones(T_fut, sys.ny);
            obj.y_max = sys.constraints.y_max * ones(T_fut, sys.ny);
            
        end
        


        function [u_fut, x_M_value, y_M_value, yalmip_error_flag] = step(obj, x0, ~, ~, y_ref)
            verbose = false; 

            %% Initialize Optimization Variables
            % Combined variables  
            u = sdpvar(obj.T_fut, obj.nu);
            y = sdpvar(obj.T_fut, obj.ny);
            x = sdpvar(obj.T_fut+1, obj.nx);

            %% Cost function
            cost = obj.Q*(y - y_ref)' * (y - y_ref) + obj.R * (u' * u);

            %% Constraints
            % Parametric model: 
            constraints = [x(1,:)' == x0];
            for k =1:1:obj.T_fut
                constraints = [constraints, x(k+1,:)' == obj.A*x(k,:)' + obj.B*u(k,:)'];
                constraints = [constraints, y(k,:)' == obj.C*x(k,:)' + obj.D*u(k,:)'];
            end

            % General control and constraints:
            constraints = [constraints;
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
            y_M_value = value(y);
            x_M_value = value(x);
        end
    end
end
