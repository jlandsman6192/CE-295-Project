%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Validation

% Param_Validation.m

clear; close all;
fs = 15;    % Font Size for plots


%% Parameter Estimates & System Matrices

%Parameter Estimates
Theta_Hat = [-0.002048298, 0.05206724, 0.017183174, 0.000184428, 0.000262303,...
    0.0772413708, 0.0026544328];
x1_eq = 69.5; %Indoor Air Temp Equilibrium [deg F]
u2_eq = 0; %Air Flow Equilibrium [CFM]
u3_eq = 71.62333; %Supply Air Temp Equilibrium [deg F]


%System matrices for identified model
Ahat = [-(Theta_Hat(1)+Theta_Hat(2)+Theta_Hat(3)+Theta_Hat(4)*u2_eq), Theta_Hat(2), Theta_Hat(3);...
    Theta_Hat(6), -(Theta_Hat(5)+Theta_Hat(6)), 0;...
    Theta_Hat(7), 0, -Theta_Hat(7)];
Bhat = [Theta_Hat(1), Theta_Hat(4)*(u3_eq-x1_eq), Theta_Hat(4)*u2_eq;...
    Theta_Hat(5), 0, 0;
    0, 0, 0];
C_dummy = eye(3);
D_dummy = 0;

% State space model
sys_hat = ss(Ahat, Bhat, C_dummy, D_dummy);

%% Load Validation Data

data = xlsread('VAV_data.xlsx');
t = data(:,1);              %time vector [hr]
air_out = data(:,2);        %outdoor air temperature, T_A [deg F]
air_supply = data(:,3);     %supply temperature, T_V [deg F]
air_in = data(:,4);         %indoor air temperature, T_Z [deg F]
mass_wall = data(:,5);      %wall mass temperature, T_W [deg F]
mass_floor = data(:,6);     %floor mass temperature, T_F [deg F]
air_flow = data(:,7);       %air flow, V [CFM]

%% Simulation

% Input vector from validation data set
U_hat = [air_out.'; air_flow.'; air_supply.'];

% Initial conditions [deg F]
That0 = [70.5; 71.9; 68.3];

% Simulate
[~,~,That] = lsim(sys_hat, U_hat, t, That0);

%% Plot Simulation

% Plot predicted and actual indoor temperature from validation data set
figure(1); clf;
plot(t, That(:,1), '-.', t, air_in)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')
