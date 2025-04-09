%% File Name: simulate_SDMPC_LTI.m
% Author: Sebastian Zieglmeier 
% Date last updated: 03.04.2025
% Description: Simulation of a LTI system controller by SD-MPC (or MPC, 
% DeePC) for a given reference trajectory
% Sources: 
% [1] - Sebastian Zieglmeier, et.al., "Semi-Data-Driven Model Pparamictive
%       Control: A Physics-Informed Data-Driven Control Approach", 
%       https://doi.org/10.48550/arXiv.2504.00746
%
%
% Notes: 
% 


close all;
clc;
clear all;


%% Get system
sys_name = "LTI_2_Tank";
sys = eval(sys_name);

%% Get scenario
scenario = 2;   % 1 - general scenario, 2 - robust scenario

%% Get controller
control_name = "rSD_MPC";   % SD_MPC, rSD_MPC, DeePC, MPC

%% Initialize Variables and control hyperparameters
T_sim = 100;       % Simulation time
T_fut = 5;         % Prediction horizon
T_d = 200;         % T_d is the number of data points used to build the Hankel matrix
T_ini = 5;         % Number of past values of DeePC

% Cost function parameters of the predictive controller
lambda_ini = 1e7;
lambda_g = 1;
r = 1e-3;
q = 1e4;


%% Load system:
% LTI state space model of the real system for the data collection
A = sys.model.A;
B = sys.model.B;
C = sys.model.C;
D = sys.model.D;

% LTI state space model of the assumed system (parametric model)
if control_name == "DeePC"
    % By nilling the parametric model the SD-MPC function can be used as
    % DeePC [1]
    A_M = zeros(size(sys.param_model.A_M));
    B_M = zeros(size(sys.param_model.B_M));
    C_M = zeros(size(sys.param_model.C_M));
    D_M = zeros(size(sys.param_model.D_M));
else
    A_M = sys.param_model.A_M;
    B_M = sys.param_model.B_M;
    C_M = sys.param_model.C_M;
    D_M = sys.param_model.D_M;
end

% LTI state space model for simulating the general or robust scenario
if scenario == 1                        % general scenario
    A_sim = sys.sim_model_1.A_sim;
    B_sim = sys.sim_model_1.B_sim;
    C_sim = sys.sim_model_1.C_sim;
    D_sim = sys.sim_model_1.D_sim;
else                                    %robust scenario
    A_sim = sys.sim_model_2.A_sim;
    B_sim = sys.sim_model_2.B_sim;
    C_sim = sys.sim_model_2.C_sim;
    D_sim = sys.sim_model_2.D_sim;
end

% System constants
nx = sys.nx; % Number of states
nu = sys.nu; % Number of inputs
ny = sys.ny; % Number of outputs

nx_M = sys.nx_M; 
nu_M = sys.nu_M; 
ny_M = sys.ny_M; 
constraints = sys.constraints;

%% Collect data for the data driven control (SD-MPC or DeePC)

u_data = zeros(T_d, nu);
u_data_sys = u_data;   
u_data_M = u_data;

y_data = zeros(T_d, ny);
y_data_sys = y_data;
y_data_M = y_data;

x_data_sys = zeros(T_d+1, nx);   
x_data_M = zeros(T_d+1, nx_M);

rand('seed', 8); % seeding for reproducibility
for i = 1:T_d
    u_data_sys(i) = rand(1)*.1*i;               % only positive values 
    [u_data_sys(i, :), w, warn] = system_boundaries(u_data_sys(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    x_data_sys(i+1, :) = (A * x_data_sys(i, :)' + B * u_data_sys(i))';
    y_data_sys(i,:) = C * x_data_sys(i, :)' + D * u_data_sys(i);
    
    u_data_M(i) = u_data_sys(i);
    y_data_M(i,:) = C_M * x_data_M(i, :)' + D_M * u_data_M(i);
    if control_name == "rSD_MPC"
        % rSD-MPC needs updating of the state of the parametric model [1]
        x_data_M(i,:) = x_data_sys(i,:);
    end
    x_data_M(i+1, :) = (A_M * x_data_M(i, :)' + B_M * u_data_M(i))';
    
    % Ensuring system constraints are satisfied:
    [x_data_sys(i+1, :), w, warn] = system_boundaries(x_data_sys(i+1, :), sys.constraints, "x");
    [y_data_sys(i, :), w, warn] = system_boundaries(y_data_sys(i, :), sys.constraints, "y");
    if w == 1
        disp(warn);
    end
end

y_data = y_data_sys - y_data_M;
u_data = u_data_M;
if control_name == "MPC"
    % Nilling the output data for using SD-MPC as pure MPC [1]
    y_data = zeros(size(y_data_sys));
end


%% Build Hankel matrizes
L = T_ini + T_fut; % Lag for right size of Hankel matrizes
num_hankel_cols = T_d - L;
H_u = u_data(1:L);
H_y = y_data(1:L);
for i = 2:num_hankel_cols+1
    H_u = [H_u, u_data(i:i+L-1)];
    H_y = [H_y, y_data(i:i+L-1)];
end





%% Intialize Controller and Simulation Variables

if control_name == "SD_MPC" || control_name == "rSD_MPC"
    control = SD_MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
elseif control_name == "MPC"
    control = MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
    % Proposition 3 in [1]:
    % control = SD_MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
elseif control_name == "DeePC"
    control = DeePC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
    % Proposition 4 in [1]:
    % control = SD_MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
else
    disp("Warning: No controller selected. Check spelling of control_name.");
end

% Matrizes to store simulation results and further values
u_sim = zeros(T_sim, 1);
x_sim = zeros(T_sim + 1, nx);
y_sim = zeros(T_sim, 1);

x_M_sim = zeros(T_sim + 1, nx_M);
y_M_sim = zeros(T_sim, 1);
x_M = zeros(T_sim + 1, nx_M);
y_M = zeros(T_sim, 1);

y_D = zeros(T_sim, 1);

%% Initial condition:
x_sim(1, :) = zeros(nx, 1);
u_past_sim = zeros(T_ini, 1);
y_past_sim = zeros(T_ini, 1);
x_0 = zeros(nx_M, 1);

%% Reference Trajectory

ref_name = "Smooth_Step_to_Sinus_paper";
ref = get_ref(ref_name, T_sim, ny, T_fut, T_ini);
% Options: Smooth_Step_to_Sinus_paper (from [1]), Smooth_Step_to_Sinus, 
% Smooth_Step, Sinus, Step
% reference trajectory values can be changed in get_ref.m
% ref = get_reference("ref_name")

%% Simulation loop:
for i = 1:T_sim 
    y_reference = ref(i:i+T_fut-1);
    [u_fut, x_M_pred, y_M_pred, opt_error_flag] = control.step(x_0, u_past_sim, y_past_sim, y_reference'); 
    y_M(i) = y_M_pred(2);
    x_M(i,:) = x_M_pred(2,:);
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim(i) = u_fut(1); % One step prediction
    % Ensure system boundaries:
    [u_sim(i, :), w, warn] = system_boundaries(u_sim(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    % Simulate real system:
    x_sim(i+1, :) = (A_sim * x_sim(i, :)' + B_sim * u_sim(i))';
    y_sim(i) = C_sim * x_sim(i, :)' + D_sim * u_sim(i);
    
    % Ensure system boundaries:
    [x_sim(i+1, :), w, warn] = system_boundaries(x_sim(i+1, :), sys.constraints, "x");
    [y_sim(i, :), w, warn] = system_boundaries(y_sim(i, :), sys.constraints, "y");
    if w == 1
        disp(warn);
    end

    % Calculate assumed model for analyse purposes
    % First calculate y_M_sim and then update state space of parametric
    % model 
    y_M_sim(i) = C_M * x_M_sim(i, :)';
    if control_name == "rSD_MPC"
        x_M_sim(i,:) = x_sim(i,:);
    end
    x_M_sim(i+1, :) = (A_M * x_M_sim(i, :)' + B_M * u_sim(i))';
    

    % Calculate data driven component 
    y_D(i) = y_sim(i) - y_M(max(1,i-1));
    % Build past data for data driven predictive control
    u_past_sim = [u_past_sim(2:end); u_sim(i)];
    y_past_sim = [y_past_sim(2:end); y_D(i)];

    % Update next state for SD-MPC as in [1]
    if control_name == "SD_MPC"
        x_0 = x_M(i,:)'; % or x_M_sim(i+1,:)'
    elseif control_name == "rSD_MPC" || control_name == "MPC"
        x_0 = x_sim(i+1,:)';
    else
        x_0 = zeros(size(x_sim(i+1,:)')); % Dummy for DeePC
    end
end


%% Numerical evaluation: 
y_ref = ref(:, 1:T_sim)';
rel_error= zeros(size(y_ref));
% Calculate relative error 
for i = 1:length(y_ref)
    if y_ref(i) ~= 0
        rel_error(i) = abs(y_sim(i) - y_ref(i)) / y_ref(i) * 100;
    else
        % If y_ref is zero, set relative error directly to zero
        rel_error(i) = 0;
    end
end
Avg_rel_error = mean(rel_error);
RMSE = (1/T_sim * sqrt((ref(:,1:T_sim)-y_sim(1:T_sim,:)')*(ref(:,1:T_sim)-y_sim(1:T_sim,:)')'));


%% Graphic evaluation: 
figure
plot(0:T_sim-1, y_sim, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Discrete timestep')
title("Output", 'FontSize', 14)
legend("y", "y_{ref}");
grid on


figure
subplot(3, 1, 1)
plot(0:T_sim-1, y_sim, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Discrete timestep')
title("Output", 'FontSize', 14)
legend("y", "y_{ref}");
grid on

subplot(3, 1, 2)
plot(0:T_sim-1, u_sim, 'g')
xlabel('Discrete timestep')
title("Control input", 'FontSize', 14)
grid on

subplot(3, 1, 3)
plot(0:T_sim-1, rel_error, 'r')
xlabel('Discrete timestep')
title("Relative Error", 'FontSize', 14)
grid on

