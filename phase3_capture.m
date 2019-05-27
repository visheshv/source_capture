%% Functions called for this script
% All except GTOC8_phase3.m, om_constants.m

%% Initialize

clear all
close all

% Mat file for the source positions 
load radiodata

% Add the path to the SNOPT directory
addpath(genpath('C:\Users\Vishesh\Documents\GTOC\GTOC 8\CD Eagle\snopt-matlab-2.5.0\'))
addpath('C:\Users\Vishesh\Documents\GTOC\GTOC 8\CD Eagle')

% define constants
global mu_earth R_earth rkcoef radiodata n_segments 
global mass_current mass_initial isp_electric

mu_earth = 398600.436233; % earth gravitational constant (DE421 value; km**3/sec**2)
R_earth  = 6378.14; % km
rkcoef =1;
n_segments= 10; % number of segments for NLP call solution
mass_current= [4e3;4e3;4e3]; % 4 mt initial mass, Not considering the change of mass
mass_initial= mass_current; % 4 mt initial mass, Not considering the change of mass
isp_electric = 5000;

% add phase change for defining true anomaly values for S2 and S3
% satellites
del_ta_2=60*pi/180;
del_ta_3=0;

%% define state vectors for S1,S2,S3 using orbit elements format
% Orbit 2 inclination of 60 deg, TA of 180 deg at Apogee
apogee_orbit_2= R_earth + 1e6; % km
perigee_orbit_2= R_earth + 384e3; % Moon position roughly 
a_orbit_2= 0.5 *(apogee_orbit_2+perigee_orbit_2); % Semi major axis [km]
e_orbit_2=(apogee_orbit_2/a_orbit_2)-1;
i_orbit_2=pi/3; 
RAAN_orbit_2=3*pi/2;
AOP_orbit_2=3*pi/2;
TA_orbit_2=pi+del_ta_2;
oev2=[a_orbit_2;e_orbit_2;i_orbit_2;AOP_orbit_2;RAAN_orbit_2;TA_orbit_2];

% Orbit 3 is the same as 2 with the RAAN flipped by 180 deg and AOP only
% pi/2
oev3 =oev2;
oev3(5)=oev2(5)-pi;
oev3(4)=pi/2;
oev3(6)=pi+del_ta_3;

% Orbit 1 is the same as Orbit2 with inclination 0 and AOP being 180 deg
% fipped
oev1 = oev2;
oev1(3)=0;
oev1(4)=pi/2;

% Convert to ECI
sv1=orb2eci(mu_earth,oev1);
sv2=orb2eci(mu_earth,oev2);
sv3=orb2eci(mu_earth,oev3);

%% Free propagation S1,S2,S3 for 3 days
% initialize initial and final times and total duration
tdur = 3; % days

sv_hist_1= propogate_earth(sv1,tdur);
sv_hist_2= propogate_earth(sv2,tdur);
sv_hist_3= propogate_earth(sv3,tdur);

[az_free,el_free]= calculate_celestial_track(sv_hist_1(:,1:3),sv_hist_2(:,1:3),sv_hist_3(:,1:3));


%% Propagation S1,S2,S3 for 3 days, check whether sources lie on track and capture
tdur = 600;
[sv_hist,capture_hist] = propagate_check(sv1,sv2, sv3,tdur);

sv_hist_1=sv_hist(:,1:3);
sv_hist_2=sv_hist(:,7:9);
sv_hist_3=sv_hist(:,13:15);

[az_1,el_1]= calculate_celestial_track(sv_hist_1,sv_hist_2,sv_hist_3);

[points,points_matrix]=points_aggregator_capture_strategy(radiodata,sv_hist_1,sv_hist_2,sv_hist_3,az_1,el_1);

% Check the points acquisition
figure(1)
plot(radiodata(:,2),radiodata(:,3),'*');
hold on
plot(az_1,el_1,'o')
hold on
plot(az_free,el_free,'r')
xlabel('Azimuth of the triangle normal (deg)')
ylabel('Elevation of the triangle normal (deg)')
title ('Free propagation vs source capturing path plan')


