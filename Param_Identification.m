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

diff_1 = zeros(10,10);
diff_2 = zeros(10,10);
diff_3 = zeros(10,10);
diff_4 = zeros(10,10);

%for i=0:9
%    for j=0:9
        % Initial conditions
        theta_hat0_1 = 10^(-5)*[1, 1, 1, 1];

        % Update Law Gain
        eta = 10^(-4);
        Gam1 = eta*eye(4);
        Gam{1} =Gam1;

        % Integrate ODEs
        [~,y] = ode23s(@(t,y) ode_gradient(t,y,data,Gam), t, theta_hat0_1);

        % Parse output
        theta_hat = y;
        
        %Difference between middle and final values
%        diff_1(i+1,j+1) = theta_hat(round(end/2),1) - theta_hat(end,1);
%        diff_2(i+1,j+1) = theta_hat(round(end/2),2) - theta_hat(end,2);
%        diff_3(i+1,j+1) = theta_hat(round(end/2),3) - theta_hat(end,3);
%        diff_4(i+1,j+1) = theta_hat(round(end/2),4) - theta_hat(end,4);
%    end
%end

%% Plot theta_hat

plot(t,theta_hat(:,1),t,theta_hat(:,2),t,theta_hat(:,3),t,theta_hat(:,4))
ylabel('theta hat','FontSize',fs)
xlabel('Time [hr]','FontSize',fs)
legend('1','2','3','4')