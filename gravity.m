function agrav = gravity(sv)

% first order equations of selenocentric orbital motion

% input

%  t  = mission elapsed time (seconds)
%  sv = state vector

% output

%  agrav = eci gravitational acceleration vector (km/sec^2)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu_earth 

r3=(norm(sv(1:3)))^3;

for i = 1:1:3
    
    agrav(i) = - mu_earth * sv(i) / r3;
    
end
