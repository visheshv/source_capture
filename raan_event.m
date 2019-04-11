function [value, isterminal, direction] = raan_event(t, y)

% twice specific orbital energy event function

% input

%  t = time (seconds)
%  y = modified equinoctial orbital elements (km, km/sec)

% output

%  value = c3 orbital energy (km/sec)^2

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global r_raan

% compute eci state vector
reci = y(1:3);
rmag = norm(reci); 

value = rmag -  r_raan;

isterminal = 1;

direction =  [];

