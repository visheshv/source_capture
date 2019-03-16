function [f, g] = capture_track(x)

% source position targeting - objective function and equality constraints

% inputs

%  x(1) = current value of transfer time (seconds)
%  x(2, nsegments) = current values of thrust histories along x,y,z

% outputs

%  f = vector of equality constraints and objective function evaluated at x

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  n_segments y_initial u_current 

% compute duration of each time interval (non-dimensional)

deltat = x(1) / n_segments;

% specify number of differential equations

neq = 18;

% truncation error tolerance

tetol = 1.0e-10;

% initialize initial time

ti = -deltat;

% set final time to current thrust duration

tof = x(1);

% set initial conditions

yi(1:18) = y_initial;

% step size guess (seconds)

h = 30.0;

% integrate for all segments

for i = 1:1:n_segments
    
    % current thrust values for Tx,Ty,Tz for S1,S2,S3
    
    u_current = x( 9 * i - 7 : 9 * i + 1 );
    
    % increment initial and final times
    
    ti = ti + deltat;
    
    tf = ti + deltat;
    
    % integrate from current ti to tf
    
    yfinal = rkf78('sel_eqm_modified', neq, ti, tf, h, tetol, yi);
    
    % reset integration vector
    
    yi = yfinal;
    
    % check for end of simulation
    
    if (tf >= tof)
        
        break;
        
    end
    
end

% objective function (thrust duration; seconds)
f(1) = x(1);

% compute inequality constraints (azimuth and elevation of the normal at the final position config)
% compute current state vector
r1 = yfinal(1:3);
r2 = yfinal(7:9);
r3 = yfinal(13:15);

[az_1,el_1]= calculate_celestial_track(r1',r2',r3');


f(2) = az_1;
f(3) = el_1;

% Thrust magnitude constraint 
Tx_S1= x( 9 * i - 7 : 9 : 9 * n_segments + 1 );
Ty_S1= x( 9 * i - 6 : 9 : 9 * n_segments + 2 );
Tz_S1= x( 9 * i - 5 : 9 : 9 * n_segments + 3 );

Tx_S2= x( 9 * i - 4 : 9 : 9 * n_segments + 4 );
Ty_S2= x( 9 * i - 3 : 9 : 9 * n_segments + 5 );
Tz_S2= x( 9 * i - 2 : 9 : 9 * n_segments + 6 );

Tx_S3= x( 9 * i - 1 : 9 : 9 * n_segments + 7 );
Ty_S3= x( 9 * i     : 9 : 9 * n_segments + 8 );
Tz_S3= x( 9 * i + 1 : 9 : 9 * n_segments + 9 );

f(4) = sum(Tx_S1.^2 + Ty_S1.^2 + Tz_S1.^2) / n_segments ;
f(5) = sum(Tx_S2.^2 + Ty_S2.^2 + Tz_S2.^2) / n_segments ;
f(6) = sum(Tx_S3.^2 + Ty_S3.^2 + Tz_S3.^2) / n_segments ;


% transpose
f = f';

% no derivatives

g = [];
