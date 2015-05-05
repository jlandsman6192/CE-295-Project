%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Identification

% Param_Identification.m

clc; clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = xlsread('VAV_data.xlsx');

%Subset data
days = 5;
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

%% Persistance of excitation

%%%% OPTION with 3-D parameter vector
phi = [(air_out-air_in), (mass_wall-air_in), (mass_floor-air_in), s]'; % <------ enter signals for 3D parametric model
t_end = t(end);
PE_mat = zeros(4);

phi_sq = zeros(4,4,length(t));
for k = 1:length(t)
    phi_sq(:,:,k) = phi(:,k) * phi(:,k)';
end
PE_mat(1,1) = 1/t_end * trapz(t, phi_sq(1,1,:));
PE_mat(2,1) = 1/t_end * trapz(t, phi_sq(2,1,:));
PE_mat(3,1) = 1/t_end * trapz(t, phi_sq(3,1,:));
PE_mat(4,1) = 1/t_end * trapz(t, phi_sq(4,1,:));
PE_mat(1,2) = 1/t_end * trapz(t, phi_sq(1,2,:));
PE_mat(2,2) = 1/t_end * trapz(t, phi_sq(2,2,:));
PE_mat(3,2) = 1/t_end * trapz(t, phi_sq(3,2,:));
PE_mat(4,2) = 1/t_end * trapz(t, phi_sq(4,2,:));
PE_mat(1,3) = 1/t_end * trapz(t, phi_sq(1,3,:));
PE_mat(2,3) = 1/t_end * trapz(t, phi_sq(2,3,:));
PE_mat(3,3) = 1/t_end * trapz(t, phi_sq(3,3,:));
PE_mat(4,3) = 1/t_end * trapz(t, phi_sq(4,3,:));
PE_mat(1,4) = 1/t_end * trapz(t, phi_sq(1,4,:));
PE_mat(2,4) = 1/t_end * trapz(t, phi_sq(2,4,:));
PE_mat(3,4) = 1/t_end * trapz(t, phi_sq(3,4,:));
PE_mat(4,4) = 1/t_end * trapz(t, phi_sq(4,4,:));

PE_lam_min = min(eig(PE_mat));
fprintf(1,'PE Level for 4D Version : %1.4f\n',PE_lam_min);


%% Problem 4(b)
% Assemble Data
data = [t, air_in, mass_wall, mass_floor, air_out, air_flow, air_supply];

% Initial conditions
multiplier = 10^(-1);
theta_hat0_1 = 2*multiplier*[1, 1, 1, 1];
theta_hat0_2 = multiplier*[2*1, 10*1];
theta_hat0_3 = 2*multiplier*[1];

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

fig1 = figure(1);
plot(t,theta_hat_1(:,1),t,theta_hat_1(:,2),t,theta_hat_1(:,3),t,theta_hat_1(:,4),...
    t,theta_hat_2(:,1),t,theta_hat_2(:,2),t,theta_hat_3(:))

ylim([-.1 1.1])
title('Progression of Parameter Identification','FontSize',fs*1.5)
xlabel('Time [hr]','FontSize',fs)
ylabel({'Value of $${\theta}$$'},'interpreter','latex','FontSize',fs)


legend({'$${\theta}_1(t)$$','$${\theta}_2(t)$$','$${\theta}_3(t)$$','$${\theta}_4(t)$$',...
    '$${\theta}_5(t)$$','$${\theta}_6(t)$$'},'interpreter','latex','FontSize',fs)

print(fig1,'.\prog_params_ID.png','-dpng');
%% Parameter Estimates & System Matrices

%Parameter Estimates
Theta_Hat = [theta_hat_1(end-2,:), theta_hat_2(end-2,:),theta_hat_3(end-2)];
% x1_eq = 69.5; %Indoor Air Temp Equilibrium [deg F]
% u2_eq = 0; %Air Flow Equilibrium [CFM]
% u3_eq = 71.62333; %Supply Air Temp Equilibrium [deg F]

Ahat = [(-Theta_Hat(1)-Theta_Hat(2)-Theta_Hat(3)), Theta_Hat(2), Theta_Hat(3);...
        Theta_Hat(6), -Theta_Hat(5)-Theta_Hat(6), 0;...
        Theta_Hat(7), 0, -Theta_Hat(7)];
    
Bhat = [Theta_Hat(1), Theta_Hat(4);...
        Theta_Hat(5), 0;...
        0, 0];

% 
% %System matrices for identified model
% Ahat = [-(Theta_Hat(1)+Theta_Hat(2)+Theta_Hat(3)+Theta_Hat(4)*u2_eq), Theta_Hat(2), Theta_Hat(3);...
%     Theta_Hat(6), -(Theta_Hat(5)+Theta_Hat(6)), 0;...
%     Theta_Hat(7), 0, -Theta_Hat(7)];
% Bhat = [Theta_Hat(1), Theta_Hat(4)*(u3_eq-x1_eq), Theta_Hat(4)*u2_eq;...
%     Theta_Hat(5), 0, 0;
%     0, 0, 0];

% Output states only (dummy variables, not used later)
C_dummy = eye(3);
D_dummy = 0;

% State space model
sys_hat = ss(Ahat, Bhat, C_dummy, D_dummy);

%% Simulation

% Input vector from validation data set
U_hat = [air_out, s];

% Initial conditions [deg F]
That0 = [air_in(1); mass_wall(1); mass_floor(1)];

% Simulate
[~,~,That] = lsim(sys_hat, U_hat, t, That0);

%% Plot Simulation

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
