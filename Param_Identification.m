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
%theta_hat0 = [0;0;0;0;0;0;0;0];

theta_hat0_1 = [0.15, 0.15, 0.15, 0.15]';
theta_hat0_2 = [0.5, 0.5]';
theta_hat0_3 = [0.5]';


theta_hat0 = [theta_hat0_1; theta_hat0_2; theta_hat0_3]';
% theta_hat0{1} = theta_hat0_1;
% theta_hat0{2} = theta_hat0_2;
% theta_hat0{3} = theta_hat0_3;

% Update Law Gain
eta = 10e-3;
Gam1 = eta*eye(4);
Gam2 = eta*eye(2);
Gam3 = eta*eye(1);

Gam{1} =Gam1;
Gam{2} =Gam2;
Gam{3} =Gam3;

% Integrate ODEs
[~,y] = ode23s(@(t,y) ode_gradient(t,y,data,Gam), t, theta_hat0);

fig1 = figure(1);
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4),t,y(:,5),t,y(:,6),t,y(:,7),'LineWidth',2); 
ylim([-.1 1.1])
title('Progression of Parameter Identification','FontSize',fs*1.5)
xlabel('Time [hr]','FontSize',fs)
ylabel({'Value of $${\theta}$$'},'interpreter','latex','FontSize',fs)


legend({'$${\theta}_1(t)$$','$${\theta}_2(t)$$','$${\theta}_3(t)$$','$${\theta}_4(t)$$',...
    '$${\theta}_5(t)$$','$${\theta}_6(t)$$'},'interpreter','latex','FontSize',fs)

print(fig1,'.\prog_params_ID.png','-dpng');

% Parse output
theta_hat = y;