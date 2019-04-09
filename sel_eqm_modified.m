function ydot = sel_eqm_modified (t, y)

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

global u_current mass_current 

% acceleration due to the moon
y1 = y(1:6);
y2 = y(7:12);
y3 = y(13:18);

u1 = u_current(1:3)/mass_current(1);
u2 = u_current(4:6)/mass_current(2);
u3 = u_current(7:9)/mass_current(3);

agrav1 = gravity(y1);
agrav2 = gravity(y2);
agrav3 = gravity(y3);
        
% compute total integration vector

ydot1 = [ y1(4)
         y1(5)
         y1(6)
         agrav1(1)+u1(1)
         agrav1(2)+u1(2)
         agrav1(3)+u1(3)];

ydot2 = [ y2(4)
         y2(5)
         y2(6)
         agrav2(1)+u2(1)
         agrav2(2)+u2(2)
         agrav2(3)+u2(3)];
ydot3 = [ y3(4)
         y3(5)
         y3(6)
         agrav3(1)+u3(1)
         agrav3(2)+u3(2)
         agrav3(3)+u3(3)];     
     
ydot = [ydot1;ydot2;ydot3];







