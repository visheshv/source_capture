clear all
close all
load radiodata

% define constants
global mu_earth R_earth rkcoef 

mu_earth = 398600.436233; % earth gravitational constant (DE421 value; km**3/sec**2)
R_earth  = 6371; % km
rkcoef =1;

% add phase change
del_ta_2=60*pi/180;
del_ta_3=0;

%% define state vectors for S1,S2,S3
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

%% Free Propagation S1 for <> year

% initialize initial and final times and total duration
tdur = 50; % days

sv_hist_1= propogate_earth(sv1,tdur);
sv_hist_2= propogate_earth(sv2,tdur);
sv_hist_3= propogate_earth(sv3,tdur);

% Plot the triangle normal on the global map
normal_1= cross(sv_hist_2(:,1:3)-sv_hist_1(:,1:3),sv_hist_3(:,1:3)-sv_hist_1(:,1:3));
normal_2= -normal_1;

[az_1, el_1, ~]= cart2sph(normal_1(:,1),normal_1(:,2),normal_1(:,3));
[az_2, el_2, ~]= cart2sph(normal_2(:,1),normal_2(:,2),normal_2(:,3));

az_1(az_1<0) = 2*pi +az_1(az_1<0);
az_2(az_2<0) = 2*pi +az_2(az_2<0);

del=0;
az_del= [(az_1*180/pi+del);(az_2*180/pi+del)];
el_del= [(el_1*180/pi+del);(el_2*180/pi+del)];

[points, points_matrix]=points_aggregator(radiodata,sv_hist_1,sv_hist_2,sv_hist_3,az_del,el_del);

% Check the points acquisition
figure(1)
plot(radiodata(:,2),radiodata(:,3),'*');
hold on
plot(az_del,el_del,'o')

figure(2)
plot3(sv_hist_1(:,1),sv_hist_1(:,2),sv_hist_1(:,3))
hold on
plot3(sv_hist_2(:,1),sv_hist_2(:,2),sv_hist_2(:,3));
hold on
plot3(sv_hist_3(:,1),sv_hist_3(:,2),sv_hist_3(:,3))





