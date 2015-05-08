%% CE 295 - Energy Systems and Control
%   Term Project
%   Non-Linear Parameter Identification
%   Prof. Moura

function xi = mle_error(Theta_Hat,t, U_hat, That0,air_in, mass_wall, mass_floor)
% This function will calculate the error between measured and simulation
% Initially set as maximum likelihood error (hence the name) but changed
% to the square error.

% This funciton is used in 'Non_linear_params_id.m'.

% Simulate building system
[That,t] = build_sim(t, U_hat, Theta_Hat, That0);

% Select parameter to optimize
y_1 = That(:,1);
y_2 = That(:,2);
y_3 = That(:,3);

% Compute the error
xi = (air_in(:)-y_1(:) + mass_wall(:)-y_2(:) + mass_floor(:)-y_3(:)).^2;

end

