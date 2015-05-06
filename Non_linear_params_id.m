%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Identification

% Param_Identification.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 10;
hours = days*24;

data = data(1:hours,:);

t = data(:,1);              %time vector [hr]
t = (0:(length(t)-1))';     %resample vector to start at 0

air_out = data(:,2);        %outdoor air temperature, T_A [deg F]
air_supply = data(:,3);     %supply temperature, T_V [deg F]
air_in = data(:,4);         %indoor air temperature, T_Z [deg F]
mass_wall = data(:,5);      %wall mass temperature, T_W [deg F]
mass_floor = data(:,6);     %floor mass temperature, T_F [deg F]
air_flow = data(:,7);       %air flow, V [CFM]

s = air_flow > 400;

figure(1)
%Plot outdoor and indoor air
plot(t,air_out,t,air_in);
ylim([50 80])
legend('Outdoor', 'Indoor')
hold on
% Plot airflow vs time
plot(t,air_flow*.2,'k');
hold off

%% Visualize inital estimates
load('./params_estimate.mat')

% Simulate inital parameters

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

% Simulation of initial estimates of parameters

% Input vector from validation data set
U_hat = [air_out, s];

% Initial conditions [deg F]
That0 = [70; 70.5; 67];

That0 = [55; 55; 55];

% Simulate
[~,~, That] = lsim(sys_hat, U_hat, t, That0);

% Plot initial parameter estimations

% Plot predicted and actual indoor temperature from validation data set
figure(1); clf;
plot(t, That(:,1), '-.', t, air_in)
title('Indoor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall)
title('Mass Wall Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor)
title('Mass Floor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

%% Parameter Estimation using lsqnonlin

%% Estimate parameters and variances
optim_options = optimset('Display', 'iter',...
'TolFun', 1e-10,... %default: 1e-4
'TolX', 1e-6,... %default: 1e-4
'MaxFunEvals', 10000,... % default: 100 * num of variables
'MaxIter', 100000,... % default: 400
'Algorithm','levenberg-marquardt'); %default: 'off'
%optim_options = [];

[p,resnorm,residual,exitflag,OUTPUT,LAMBDA,Jacobian] = lsqnonlin(@mle_error, Theta_Hat, [],[],optim_options,...
    t, U_hat, That0, air_in, mass_wall, mass_floor);
disp(' ')
p

save('p_values.mat','p');

load('Best_p_values.mat');

Theta_Hat = p;
% Simulate inital parameters

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

% Simulation of initial estimates of parameters

% Input vector from validation data set
U_hat = [air_out, s];

% Initial conditions [deg F]
That0 = [70; 70.5; 67];

%That0 = [55; 55; 55];

% Simulate
[~,~, That] = lsim(sys_hat, U_hat, t, That0);

% Plot initial parameter estimations

% Plot predicted and actual indoor temperature from validation data set
fig1 = figure(1); clf;
plot(t, That(:,1), '-.', t, air_in)
title('Indoor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

print(fig1,'.\nonlin_estimate_air_in.png','-dpng');

fig2 = figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall)
title('Mass Wall Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Save plot 
print(fig2,'.\nonlin_estimate_mass_wall.png','-dpng');

fig3 = figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor)
title('Mass Floor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Save plot 
print(fig3,'.\nonlin_estimate_mass_floor.png','-dpng');

