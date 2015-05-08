%% CE 295 - Energy Systems and Control
%   Term Project
%   Optimization for states
%   Prof. Moura


% Optimization.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 10;
hours = days*24;

% Times for training data
t_0 = [10:hours];    % training data

% Times for validation data
t_1 = [400:600];    % validation data 1; state is always 0
t_2 = [1575:1750];  % validation data 2; state has night ventilation 

data = data(t_2,:);

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

% Comfort Range
T_cof_lb = zeros(length(occ),1);
T_cof_ub = zeros(length(occ),1);

T_cof_lb(occ) = 69;
T_cof_lb(~occ) = 60; 

T_cof_ub(occ) = 78;
T_cof_ub(~occ) = 85; 

% Figure out the different states from air_flow
s = air_flow > 400;


figure(1)
%Plot outdoor and indoor air
plot(t,air_out,t,air_in,'LineWidth',1.5);
ylim([50 80])
title('Building Trends','FontSize',fs*1.5)
ylabel('Temperature','FontSize',fs);
xlabel('Time [hr]','FontSize',fs);
hold on
% Plot airflow vs time
plot(t,air_flow*.2,'k','LineWidth',1.5);
legend('Outdoor', 'Indoor','Airflow');
hold off

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

% Plot simulation trends

% Plot initial parameter estimations
T_op = (That(:,1) + That(:,2) + That(:,3))./3;
T_opact = (air_in + mass_wall + mass_floor)./3;

fig = figure(1); clf;
plot(t, T_op, '-.','LineWidth',2.5)         % Predicted Operative Temp
hold on
plot(t, T_opact, 'b','LineWidth',1.5)        % Actual Operative Temp
plot(t, T_cof_lb,'--g','LineWidth',1.5);    % lower bound
plot(t, T_cof_ub,'--g','LineWidth',1.5);    % upper bound

plot(t, 60*s,'r','LineWidth',1.5);          % state variable;

plot(t,air_out,'k')                         % Outdoor air

ylim([55 90]);
hold off

title('Operative Temperature Prediction','FontSize',fs*1.5)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted Operative Temp','True Operative Temp','','Comfort Bounds','State Variable','Outside Air Temp')

% Save Plot
%print(fig,'.\test_2_Top_.png','-dpng');

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

s_subset = 17:32;
S = s(s_subset);
air_out_sub = air_out(17:32);
t_s = t(s_subset);

%Initialize simulation

U_hat = [air_out_sub, S];       % Input vector from validation data set
That0 = [70; 70.5; 67];     % Initial conditions [deg F]

% Simulate
[~,~, That] = lsim(sys_hat, U_hat, t_s, That0);

% Plot initial parameter estimations
T_op = (That(:,1) + That(:,2) + That(:,3))./3;

% Add offset to predicted T_op to continue with optimization
T_op = T_op + 5;

fig = figure(1); clf;
plot(t_s, T_op, '-.','LineWidth',2.5)                   % Predicted Operative Temp
hold on
plot(t_s, T_cof_lb(s_subset),'--g','LineWidth',1.5);    % lower bound
plot(t_s, T_cof_ub(s_subset),'--g','LineWidth',1.5);    % upper bound

plot(t_s, 60*S,'r','LineWidth',1.5);                    % state variable;

plot(t_s,air_out,'k')                                   % Outdoor air

ylim([55 90]);
hold off

title('Operative Temperature Prediction','FontSize',fs*1.5)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted Operative Temp','Comfort Bounds','','State Variable','Outside Air Temp')

