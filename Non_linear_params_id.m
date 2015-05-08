%% CE 295 - Energy Systems and Control
%   Term Project
%   Non-Linear Parameter Identification
%   Prof. Moura
%
%   This script will perform a non linear optimization of parameters.
%   Initial parameters are taken from the final results from gradient 
%   descent method.  

%   Functions need to run this script are 'mle_error.m' and 'build_sim.m'.

% Non_linear_param_id.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 10;
hours = days*24;

% Times for training data
t_0 = [10:hours];           % training data

% Times for validation data
t_1 = [400:600];            % validation data 1; state is always 0
t_2 = [1575:1750];          % validation data 2; state has night ventilation 

data = data(t_0,:);         % Subset data

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

% %Convert problem into a four state problem
% s_2 = air_flow <= 400 & air_flow > 200;
% s_3 = air_flow <= 200;
% 
% s_fin = air_flow;
% 
% s_fin(s_3)=0;
% s_fin(s_2)=1;
% s_fin(s)=4;
% 
% s = s_fin;

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

% Plot state variable of the data
fig2 = figure(2);
plot(t,s,'LineWidth',1.5);
ylim([0 2]);
title('State variable','FontSize',fs*1.5)
ylabel('State','FontSize',fs);
xlabel('Time [hr]','FontSize',fs);

%print(fig2,'.\state_variable.png','-dpng');

%% Visualize inital parameter estimates
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

% Simulate
[~,~, That] = lsim(sys_hat, U_hat, t, That0);

% Plot initial parameter estimations

% Plot indoor air temp predicted results with actual results
figure(1); clf;
plot(t, That(:,1), '-.', t, air_in,'LineWidth',1.5)
title('Indoor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Plot mass wall temp predicted results with actual results
figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall,'LineWidth',1.5)
title('Mass Wall Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Plot mass floor temp predicted results with actual results
figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor,'LineWidth',1.5)
title('Mass Floor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

%% Parameter Estimation using lsqnonlin

% Options for for the algorithm
optim_options = optimset('Display', 'iter',...
'TolFun', 1e-10,... %default: 1e-4
'TolX', 1e-6... %default: 1e-4
);
%'MaxFunEvals', 700,... % default: 100 * num of variables
%'MaxIter', 100000,... % default: 400
%'Algorithm','levenberg-marquardt'); %default: 
%optim_options = [];

% Nonlinear optimization of parameters
[p,resnorm,residual] = lsqnonlin(@mle_error, Theta_Hat, [],[],optim_options,...
    t, U_hat, That0, air_in, mass_wall, mass_floor);
disp('Finished parameter estimation and new theta hat is')
p

% Save p values for later use
% uncomment when performing predictions
%save('p_values.mat','p');

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

% Plot final parameter estimations

% Plot predicted and actual indoor temperature from validation data set
fig1 = figure(1); clf;
plot(t, That(:,1), '-.', t, air_in,'LineWidth',1.5)
title('Indoor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Save Plot
%print(fig1,'.\nonlin_estimate_air_in.png','-dpng');

% Plot predicted and actual mass wall temperature from validation data set
fig2 = figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall,'LineWidth',1.5)
title('Mass Wall Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Save plot 
%print(fig2,'.\nonlin_estimate_mass_wall.png','-dpng');

% Plot predicted and actual mass floortemperature from validation data set
fig3 = figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor,'LineWidth',1.5)
title('Mass Floor Temperature Prediction','FontSize',fs*1.5)
ylim([55 75]);
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

% Save plot 
%print(fig3,'.\nonlin_estimate_mass_floor.png','-dpng');

