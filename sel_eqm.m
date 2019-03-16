function ydot = sel_eqm (t, y)

% selenocentric equations of motion

% includes perturbations due to:

%   non-spherical lunar gravity
%   sun and earth point-mass gravity

% input

%  t = simulation time (seconds)
%  y = state vector (kilometers and kilometers/second)

% output

%  ydot = integration vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% acceleration due to the moon

agrav = gravity(y);
        
% compute total integration vector

ydot = [ y(4)
         y(5)
         y(6)
         agrav(1)
         agrav(2)
         agrav(3)];








