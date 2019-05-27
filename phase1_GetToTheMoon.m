%% Get to the Moon (Phase 1) 
% Reference: GTOC8 Results and Methods of ESA Advanced Concepts Team and
% JAXA-ISAS

clear all 
close all

global T_engine Isp_engine emu r_raan rkcoef nu_inv_control 
global r0 v0 m0 omega_sat dt_RAAN

rkcoef=1;

%% Define initial states of the satellites
t0 = 58849.0; % MJD
R_earth  = 6378.14; % km
emu = 398600.436233; % kg m3/ s2
r0 = [-R_earth+400; 0; 0]; % Initial position km, for S1,S2,S3
v0 = [0; -sqrt(emu/norm(r0)); 0]; % Initial velocity km/s for S1,S2,S3
m0 = 4e3; 
omega_sat = norm(v0)/norm(r0); % rad/s

%% Define propulsion system characteristics
T_engine   = 0.1; 
Isp_engine = 5000;

%% Define moon position and epoch
% epoch  time 
t0 = 58849.0; % MJD
a= 383500.0; % semimajor axis, a (km) 
e= 0.04986;  % eccentricity, e 
i= 5.2586;   % inclination, i (deg.) 
RAAN = 98.0954; % LAN, (deg.) 
AOP = 69.3903; % Arg. peri., ! (deg.) 
M0= 164.35025; %Mean anomaly, M0 (deg.) 
mu_moon = 4902.8006; % kg m3 / s2

%% Calculate position magnitude of the ascending node
nu_raan = 360 - AOP;
rp_moon = a *(1-e);
h = sqrt(mu_moon * rp_moon * (1+e)); % Source: http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node33.html
r_raan = (h^2 / mu_moon ) / (1 + e* cosd(nu_raan) ); % Source: https://www.physicsforums.com/threads/position-of-ascending-and-descending-nodes.534665/

% find time to RAAN
E= mod(2 * atan3(tand(nu_raan * 0.5) , sqrt(1+e) /sqrt(1-e)),2*pi);
M_raan=E-e*sin(E);
dt_RAAN= (M_raan-M0*pi/180) * sqrt(a^3/emu);


%% Integrate the position of the satellite until it reaches |r|>r_raan adn find root such that there exists a rendezvous with the moon at the ascending node

nu_inv_control = 60 * pi/180;

T0 = [0.8*T_engine T_engine];
func = @residual_calc_phase1;
T_root = fzero(func,T0);

[rf, tf, capture_states] = inverse_control_law_integration(r0,v0,m0, T_root,nu_inv_control);

%Plot the states
figure(1)
plot3(capture_states(:,1),capture_states(:,2), capture_states(:,3));
xlabel('x position (km)');
ylabel('y position (km)');
zlabel('z position (km)');

figure(2);
plot((capture_states(:,4).^2+capture_states(:,5).^2+ capture_states(:,6).^2).^0.5);
ylabel('Speed (km/s)')

figure(3);
plot(capture_states(:,7));
ylabel('Mass (kg)')
