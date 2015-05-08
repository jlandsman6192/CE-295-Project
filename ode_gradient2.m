%% CE 295 - Energy Systems and Control
%   Term Project
%   Parameter Identification
%   Prof. Moura

% ode_gradient2.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% theta_h   : parameter estimate
% data      : input-output data used to feed algorithm
% Gam       : Update law gain

function [theta_h_dot2] = ode_gradient2(t,theta_h,data,Gam)


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

% Compute mass wall temperature at NEXT time step
mass_wall_plus = interp1(it,imass_wall,t+dt);

% Compute \dot{T} using forward difference in time 
% z = \dot{T} = (T(t+dt) - T(t))/dt
z = (mass_wall_plus - mass_wall)/dt;

% Assemble regressor vector, \phi
phi = [(air_out-mass_wall), (air_in-mass_wall)]';

%% Gradient Update Law
% Normalization signal
msq = 1 + phi'*phi;

% Estimation error: \epsilon = z - \theta_h^T \phi
epsilon = (z - theta_h'*phi)/msq;

% Update Law
theta_h_dot2 = Gam * epsilon * phi;