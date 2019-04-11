function [rf, tf,capture_states] = inverse_control_law_integration(r0,v0,m0, T,nu_inv_control)

global T_control nu_inversion_control r_raan

T_control  = T;
nu_inversion_control = nu_inv_control;

tof = 20 * 86400; % Assume 20 days for intersection
dt = 0.001 * 86400;
yi = [r0;v0;m0];

% number of equations
neq = 7;
% step size guess (seconds)
h = 30;
tetol = 1.0e-3;
ti = -dt;
i=0;

while(1)
    
    % increment counter
    i=i+1;
    % increment initial and final times
    ti = ti + dt;
    tf = ti + dt;
    
    % integrate from current ti to tf
    yfinal = rkf78('eqm_ICL_phase1', neq, ti, tf, h, tetol, yi);
    yi = yfinal;
    capture_states(i,:) = yi';
    
    r_temp =norm(yfinal(1:3));
    % check for end of simulation
    if (tf >= tof) || r_temp> r_raan
        
        break;
        
    end
end

rf = yfinal(1:3);
end
