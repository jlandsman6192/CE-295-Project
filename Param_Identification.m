%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Identification

% Param_Identification.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');


t = data(:,1);              %time vector [hr]
air_out = data(:,2);        %outdoor air temperature, T_A [deg F]
air_supply = data(:,3);     %supply temperature, T_V [deg F]
air_in = data(:,4);         %indoor air temperature, T_Z [deg F]
mass_wall = data(:,5);      %wall mass temperature, T_W [deg F]
mass_floor = data(:,6);     %floor mass temperature, T_F [deg F]
air_flow = data(:,7);       %air flow, V [CFM]

%% Problem 4(b)
% Assemble Data
data = [t, air_in, mass_wall, mass_floor, air_out, air_flow, air_supply];

% Initial conditions
multiplier = 10^(-1);
theta_hat0_1 = multiplier*[1, 1, 1, 1];
theta_hat0_2 = multiplier*[1, 1];
theta_hat0_3 = multiplier*[1];

% Update Law Gain
eta = 10^(-2);
Gam1 = eta*eye(4);
Gam2 = eta*eye(2);
Gam3 = eta*eye(1);

% Integrate ODEs
[~,y1] = ode23s(@(t,y) ode_gradient1(t,y,data,Gam1), t, theta_hat0_1);
[~,y2] = ode23s(@(t,y) ode_gradient2(t,y,data,Gam2), t, theta_hat0_2);
[~,y3] = ode23s(@(t,y) ode_gradient3(t,y,data,Gam3), t, theta_hat0_3);

% Parse output
theta_hat_1 = y1;
theta_hat_2 = y2;
theta_hat_3 = y3;        

%% Plot theta_hat

plot(t,theta_hat_1(:,1),t,theta_hat_1(:,2),t,theta_hat_1(:,3),t,theta_hat_1(:,4),...
    t,theta_hat_2(:,1),t,theta_hat_2(:,2),t,theta_hat_3(:))
ylabel('theta hat','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('1','2','3','4','5','6','7')