function [That,t] = build_sim(t, U_hat, Theta_Hat, That0)
% This function will simulation the results of the system

Ahat = [(-Theta_Hat(1)-Theta_Hat(2)-Theta_Hat(3)), Theta_Hat(2), Theta_Hat(3);...
        Theta_Hat(6), -Theta_Hat(5)-Theta_Hat(6), 0;...
        Theta_Hat(7), 0, -Theta_Hat(7)];
    
Bhat = [Theta_Hat(1), Theta_Hat(4);...
        Theta_Hat(5), 0;...
        0, 0];
    
% Output states only (dummy variables, not used later)
C_dummy = eye(3);
D_dummy = 0;

sys_hat = ss(Ahat, Bhat, C_dummy, D_dummy);
    
% Simulate
[That,t] = lsim(sys_hat, U_hat, t, That0);

end

