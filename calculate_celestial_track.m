function [az,el] = calculate_celestial_track(sv_hist_1,sv_hist_2,sv_hist_3)

% Calculate normal celestial track

% Plot the triangle normal on the global map
normal_1= cross(sv_hist_2(:,1:3)-sv_hist_1(:,1:3),sv_hist_3(:,1:3)-sv_hist_1(:,1:3));

[az, el, ~]= cart2sph(normal_1(:,1),normal_1(:,2),normal_1(:,3));

% Rad to deg conversion
az(az<0) = 2*pi +az(az<0);
az = (180/pi) * az;  % Achieved source position
el = (180/pi) * el;