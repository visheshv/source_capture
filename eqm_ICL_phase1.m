function ydot = eqm_ICL_phase1 (t, y)

% Equations of motion

% Inversion control law

% input

%  t = simulation time (seconds)
%  y = state vector (kilometers and kilometers/second and kg)

% output

%  ydot = integration vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global T_control Isp_engine emu nu_inversion_control

% acceleration due to the moon
r = y(1:3);
v = y(4:6);
m = y(7);

oev = eci2orb1 (emu, r, v);

f = oev(6); % True anomaly

gravity = - emu * r / norm(r)^3;

% Fixed perigee control law
if f > (pi - nu_inversion_control)  && f < (pi + nu_inversion_control) 
    u = -v / norm(v); 
else
    u = v / norm(v);
end

rdot = v;
vdot = gravity + (T_control / m) * u;
mdot = -T_control / (Isp_engine * 9.806);

ydot = [rdot;vdot;mdot];

