%% CE 295 - Energy Systems and Control
%   HW 2 : Parameter Identification for a Smart Home Thermostat
%   Jared Landsman, SID 25952611
%   Prof. Moura

% ode_gradient.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% theta_h   : parameter estimate
% data      : input-output data used to feed algorithm
% Gam       : Update law gain

function [theta_h_dot] = ode_gradient(t,theta_h,data,Gam)


%% Parse Input Data
it = data(:,1);             %time vector [hr]
iair_in = data(:,2);        %outdoor air temperature [deg F]
imass_wall = data(:,3);     %supply temperature [deg F]
imass_floor = data(:,4);    %indoor air temperature [deg F]
iair_out = data(:,5);       %wall mass temperature [deg F]
iair_flow = data(:,6);      %floor mass temperature [deg F]
iair_supply = data(:,7);    %air flow [CFM]

%% Interpolate data
air_in = interp1(it,iair_in,t);
mass_wall = interp1(it,imass_wall,t);
mass_floor = interp1(it,imass_floor,t);
air_out = interp1(it,iair_out,t);
air_flow = interp1(it,iair_flow,t);
air_supply = interp1(it,iair_supply,t);

%% Parametric model notation
% Samping time step
dt = 1;

% Compute Room temperature at NEXT time step
air_in_plus = interp1(it,iair_in,t+dt);
%mass_wall_plus = interp1(it,imass_wall,t+dt);
%mass_floor_plus = interp1(it,imass_floor,t+dt);

% Compute \dot{T} using forward difference in time 
% z = \dot{T} = (T(t+dt) - T(t))/dt
z1 = (air_in_plus - air_in)/dt;
%z2 = (mass_wall_plus - mass_wall)/dt;
%z3 = (mass_floor_plus - mass_floor)/dt;

% Assemble regressor vector, \phi
%phi = [air_in; mass_wall; mass_floor; air_out; air_flow; air_supply];

phi_1 = [(air_out-air_in), (mass_wall-air_in), (mass_floor-air_in), air_flow.*(air_supply-air_in)*60]';
%phi_2 = [(air_out-mass_wall), (air_in-mass_wall)]';
%phi_3 = [(air_in-mass_floor)]';

%% Gradient Update Law
% Normalization signal
msq_1 = 1 + phi_1'*phi_1;
%msq_2 = 1 + phi_2'*phi_2;
%msq_3 = 1 + phi_3'*phi_3;


% Estimation error: \epsilon = z - \theta_h^T \phi
epsilon1 = (z1 - theta_h(1:4)'*phi_1)/msq_1;
%epsilon2 = (z2 - theta_h(5:6)'*phi_2)/msq_2;
%epsilon3 = (z3 - theta_h(7)'*phi_3)/msq_3;

% Update Law
theta_h_dot1 = Gam{1} * epsilon1 * phi_1;
%theta_h_dot2 = Gam{2} * epsilon2 * phi_2;
%theta_h_dot3 = Gam{3} * epsilon3 * phi_3;

% theta_h_dot{1} = theta_h_dot1;
% theta_h_dot{2} = theta_h_dot2;
% theta_h_dot{3} = theta_h_dot3;

%theta_h_dot = [theta_h_dot1; theta_h_dot2; theta_h_dot3];

theta_h_dot = theta_h_dot1;