function theta_residual = residual_calc_phase1(T)

global nu_inv_control r0 v0 m0 omega_sat dt_RAAN

[rf, tf, capture_states] = inverse_control_law_integration(r0,v0,m0, T,nu_inv_control);
theta2 = atan3(rf(2),rf(1));
theta1 = mod(omega_sat * (tf - dt_RAAN),2*pi);

theta_residual = theta1-theta2;