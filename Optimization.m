%% CE 295 - Energy Systems and Control
%   Term Project
%   Non-Linear Parameter Identification

% Non_linear_param_id.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 10;
hours = days*24;

% Times for validation data
t_0 = [10:hours];    % training data
t_1 = [400:600];    % validation data 1; state is always 0
t_2 = [1575:1750];  % validation data 2; state has night ventilation 

data = data(t_0,:);

t = data(:,1);              %time vector [hr]
t = (0:(length(t)-1))';     %resample vector to start at 0

air_out = data(:,2);        %outdoor air temperature, T_A [deg F]
air_supply = data(:,3);     %supply temperature, T_V [deg F]
air_in = data(:,4);         %indoor air temperature, T_Z [deg F]
mass_wall = data(:,5);      %wall mass temperature, T_W [deg F]
mass_floor = data(:,6);     %floor mass temperature, T_F [deg F]
air_flow = data(:,7);       %air flow, V [CFM]
hour = data(:,8);           %time of day in [HH]

% Decide whether time is unoccupied or occupied
occ = hour >= 8 & hour <= 17;

% Figure out the different states from air_flow
s = air_flow > 400;


%% Parameter Estimation using lsqnonlin

% Load best values for parameters
load('Best_p_values.mat');

Theta_Hat = p;

% Simulate final parameters
Ahat = [(-Theta_Hat(1)-Theta_Hat(2)-Theta_Hat(3)), Theta_Hat(2), Theta_Hat(3);...
        Theta_Hat(6), -Theta_Hat(5)-Theta_Hat(6), 0;...
        Theta_Hat(7), 0, -Theta_Hat(7)];
    
Bhat = [Theta_Hat(1), Theta_Hat(4);...
        Theta_Hat(5), 0;...
        0, 0];
    
% Output states only (dummy variables, not used later)
C_dummy = eye(3);
D_dummy = 0;

% State space model
sys_hat = ss(Ahat, Bhat, C_dummy, D_dummy);

% Simulation of final estimates of parameters

% Input vector from validation data set
U_hat = [air_out, s];

% Initial conditions [deg F]
That0 = [70; 70.5; 67];

% Simulate
[~,~, That] = lsim(sys_hat, U_hat, t, That0);

%% Optimization of the state variable

% Grid State and Preallocate
s_grid = [0 1]';

% Grid size
ns = 2;  % No. of states

% Planning horizon (time steps)
N = length(t);

% Preallocate Value Function (rows index state, columns index time)
V = inf*ones(ns,N+1);

% Preallocate Control (rows index state, columns index time)
u_star = zeros(ns,N);

%% Solve DP Version 1
tic;

% Boundary Condition of Value Function (Principle of Optimality)
V(:,N+1) = 0;

% Iterate time steps
for k = N:-1:1

    % Iterate over grid
    for idx = 1:ns
        
        % Operative Temp dominant bounds
        lb = 60; ub = 85;
        T_op_grid = linspace(lb,ub,50)';
        
        % Calculate next s (vectorized for all T_op_grid)
        for i = 1:length(T_op_grid)
            if k > 1
                if T_op_grid(i)<=69 && T_op_grid(i)>=78 && occ(k-1)==1
                    s_nxt(i) = 1;
                else
                    s_nxt(i) = 0;
                end
            else
                s_nxt(i) = 0;
            end
        end
        
        % Compute value function at nxt time step (need to interpolate)
        V_nxt = interp1(s_grid,V(:,k+1),s_nxt,'linear');
        
        % Value Function (Principle of Optimality)
        [V(idx, k), ind] = min(V_nxt);
        
        % Save Optimal Control
        u_star(idx,k) = T_op_grid(ind);

    end
end

solveTime = toc;
fprintf(1,'DP Solver Time %2.2f sec \n',solveTime);

%% Optimization Method 2

N = length(t);

% Compute all possible operative temperatures
for k = 1:N
       
        % Initial conditions
        if k == 1
            That0_1 = [70; 70.5; 67];
            That0_2 = [70; 70.5; 67];
            That0_3 = [70; 70.5; 67];
            That0_4 = [70; 70.5; 67];
        else
            That0_1 = [That_air(k-1,1);That_wall(k-1,1);That_floor(k-1,1)];
            That0_2 = [That_air(k-1,2);That_wall(k-1,2);That_floor(k-1,2)];
            That0_3 = [That_air(k-1,3);That_wall(k-1,3);That_floor(k-1,3)];
            That0_4 = [That_air(k-1,4);That_wall(k-1,4);That_floor(k-1,4)];
        end
        
        %Inputs
        U_hat_1 = [air_out(k-1:k), [0;0] ];
        U_hat_2 = [air_out(k-1:k), [0;1] ];
        U_hat_3 = [air_out(k-1:k), [1;0] ];
        U_hat_4 = [air_out(k-1:k), [1;1] ];
        
        %Simulations
        [~,~, That_opt_11] = lsim(sys_hat, U_hat_1, t(k-1:k), That0_1);
        [~,~, That_opt_13] = lsim(sys_hat, U_hat_1, t(k-1:k), That0_3);
        [~,~, That_opt_21] = lsim(sys_hat, U_hat_2, t(k-1:k), That0_1);
        [~,~, That_opt_23] = lsim(sys_hat, U_hat_2, t(k-1:k), That0_3);
        [~,~, That_opt_32] = lsim(sys_hat, U_hat_3, t(k-1:k), That0_2);
        [~,~, That_opt_34] = lsim(sys_hat, U_hat_3, t(k-1:k), That0_4);
        [~,~, That_opt_42] = lsim(sys_hat, U_hat_4, t(k-1:k), That0_2);
        [~,~, That_opt_44] = lsim(sys_hat, U_hat_4, t(k-1:k), That0_4);
        
        %Collect temperatures
        That_air(k,1:4) = [That_opt_1(2,1) That_opt_2(2,1) That_opt_3(2,1) That_opt_4(2,1)];
        That_wall(k,1:4) = [That_opt_1(2,2) That_opt_2(2,2) That_opt_3(2,2) That_opt_4(2,2)];
        That_floor(k,1:4) = [That_opt_1(2,3) That_opt_2(2,3) That_opt_3(2,3) That_opt_4(2,3)];
        
        %Calculate operative temperature
        for i = 1:4
            T_op(k,i) = 1/3*(That_air(k,i) + That_wall(k,i) + That_floor(k,i));
        end
end