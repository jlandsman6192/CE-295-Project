%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Identification

% Param_Identification.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 14;
hours = days*24;
starting_hour = 1;

data = data(starting_hour:hours,:);

t = data(:,1);              %time vector [hr]
t = (0:(length(t)-1))';     %resample vector to start at 0

air_out = data(:,2);        %outdoor air temperature, T_A [deg F]
air_supply = data(:,3);     %supply temperature, T_V [deg F]
air_in = data(:,4);         %indoor air temperature, T_Z [deg F]
mass_wall = data(:,5);      %wall mass temperature, T_W [deg F]
mass_floor = data(:,6);     %floor mass temperature, T_F [deg F]
air_flow = data(:,7);       %air flow, V [CFM]

s = air_flow > 400;

%% Parameter Estimates & System Matrices

%Room properties
Room_L = 38; %ft
Room_W = 49.5; %ft
Room_H = 10; %ft ESTIMATE
Thick_Conc = 4; %in
Thick_Cem = 2; %in
Thick_Ceil = .5; %in ESTIMATE
[max_air_flow,idx] = max(air_flow); %CFM

%Material properties
rho_Conc = 145; %lb/ft^3
rho_Cem = 95; %lb/ft^3
rho_Air = 32.174*2.329e-3; %lb/ft^3
Cp_Conc = .23; %BTU/(lb*deg F)
Cp_Cem = .37; %BTU/(lb*deg F)
Cp_Air = .2403; %BTU/(lb*deg F)
R_Conc = .07; %deg F*ft^2*hr/(BTU*in)
R_Cem = .26; %deg F*ft^2*hr/(BTU*in)
R_Ceil = .45; %deg F*ft^2*hr/(BTU*in) ESTIMATE
R_Ins = 20; %deg F*ft^2*hr/BTU
Film_In = .68; %deg F*ft^2*hr/BTU
Film_Out = .17; %deg F*ft^2*hr/BTU

%Real parameters
P = rho_Air*Cp_Air*max_air_flow*60*(air_out(idx)-air_in(idx)); %BTU/hr
R_AZ = (R_Ceil*Thick_Ceil + Film_In + Film_Out)/(Room_L*Room_W); %deg F * hr/BTU
R_FZ = (R_Conc*Thick_Conc)/(Room_L*Room_W); %deg F * hr/BTU
R_WZ = (R_Cem*Thick_Cem + Film_In)/((Room_L + Room_W)*Room_H*2); %deg F * hr/BTU
R_AW = (R_Ins + Film_Out)/((Room_L + Room_W)*Room_H); %deg F * hr/BTU
C_Z = rho_Air*Cp_Air*(Room_L*Room_W*Room_H); %BTU/deg F
C_W = rho_Cem*Cp_Cem*(Room_L + Room_W)*Room_H*2*Thick_Cem/12; %BTU/deg F
C_F = rho_Conc*Cp_Conc*(Room_L*Room_W*Thick_Conc/12); %BTU/deg F

real_params = [1/(C_Z*R_AZ) 1/(C_Z*R_WZ) 1/(C_Z*R_FZ) P...
    1/(C_W*R_AW) 1/(C_W*R_WZ) 1/(C_F*R_FZ)];

%System matrices for real parameters
Ahat = [-1/C_Z*(1/R_AZ + 1/R_WZ + 1/R_FZ) 1/(C_Z*R_WZ) 1/(C_Z*R_FZ);...
    1/(C_W*R_WZ) -1/C_W*(1/R_AW + 1/R_WZ) 0;...
    1/(C_F*R_FZ) 0 -1/(C_F*R_FZ)];
Bhat = [1/(C_Z*R_AZ) P/C_Z;...
    1/(C_W*R_AW) 0; 0 0];
C_dummy = eye(3);
D_dummy = 0;

% State space model
sys_hat = ss(Ahat, Bhat, C_dummy, D_dummy);

%% Simulation

% Input vector from validation data set
U_hat = [air_out.'; s.'];

% Initial conditions [deg F]
That0 = [air_in(1); mass_wall(1); mass_floor(1)];

% Simulate
[~,~,That] = lsim(sys_hat, U_hat, t, That0);

%% Plot Simulation

% Plot predicted and actual indoor temperature from validation data set
figure(1); clf;
plot(t, That(:,1), '-.', t, air_in,'LineWidth',1.5)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(2); clf;
plot(t, That(:,2), '-.', t, mass_wall,'LineWidth',1.5)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')

figure(3); clf;
plot(t, That(:,3), '-.', t, mass_floor,'LineWidth',1.5)
ylabel('Temperature [deg F]','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('Predicted','True')
