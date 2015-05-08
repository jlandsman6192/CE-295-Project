%% CE 295 - Energy Systems and Control
%   HW 2 : Parameter Identification for a Smart Home Thermostat
%   Jared Landsman, SID 25952611
%   Prof. Moura

% ode_gradient1.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% theta_h   : parameter estimate
% data      : input-output data used to feed algorithm
% Gam       : Update law gain

function [theta_h_dot1] = ode_gradient1(t,theta_h,data,Gam)


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

% Figure out the different states from air_flow
s = air_flow > 400;

%% Parametric model notation
% Samping time step
dt = 1;

% Compute Room temperature at NEXT time step
air_in_plus = interp1(it,iair_in,t+dt);

% Compute \dot{T} using forward difference in time 
% z = \dot{T} = (T(t+dt) - T(t))/dt
z = (air_in_plus - air_in)/dt;

% Assemble regressor vector, \phi
phi = [(air_out-air_in), (mass_wall-air_in), (mass_floor-air_in), s]';

%% Gradient Update Law
% Normalization signal
msq = 1 + phi'*phi;

% Estimation error: \epsilon = z - \theta_h^T \phi
epsilon = (z - theta_h'*phi)/msq;

% Update Law
theta_h_dot1 = Gam * epsilon * phi;