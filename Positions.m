close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D'); % Working directory 

%% PARAMETERS

mirr2view_dist = 0.585; % Distance between the mirror and the viewport [m]
view_thick = 0.06; % Viewport thickness [m]
n_air = 1.000293; % Air refractive index
n_viewport = 1.46; % Fused silica refractive index 
n_gas = 1.000293; % Gas mixture refractive index, here taken as air [To Be Modified]

sphere_rad = 0.10; % Sphere radius 
Am = 0.0075; % Curvature radius [m]

angle = zeros(1200,1); % Intrapolated mirror positions based on the recorded data

%% POSITION DEVIATION COMPUTATION 

% struct_runs = dir([WorkDir '\*.txt']);
% nb_runs = length(struct_runs);
% runs_list = sortrows(char(struct_runs.name));
% 
% positions = zeros(1200,9);
% 
% for i = 1:nb_runs
%     data = readtable([WorkDir '\' runs_list(i,:)]);
%     positions(:,i) = data.Output_um_;
% end
% 
% time = data.time_ms_;
% position_mean = mean(positions,2); % mean mirror position [mrad]
% position_std = std(positions,0,2);
% 
% figure(1)
% hold on
% fill([time', fliplr(time')], [(position_mean+position_std)', fliplr((position_mean-position_std)')], 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
% plot(time, position_mean, 'k-')
% hold off
% xlabel('Time [ms]')
% ylabel('Mirror position [urad]')

data_mirror = readtable([WorkDir '\Positions.txt']);
mirror_time = data_mirror.time_ms_; % Position sampling time [ms]
mirror_angle = data_mirror.Output_um_; % Mirror angle [mrad]

%% MIRROR POSITION APPROXIMATION

time = 0.025:50e-3:61; % Reconstructed time vector [ms]

for i = 1:length(time)
    tminus_index = find(mirror_time<time(i),1,'last');
    tplus_index = find(mirror_time>time(i),1,'first');
    tminus = mirror_time(tminus_index);
    tplus = mirror_time(tplus_index);
    anglminus = mirror_angle(tminus_index);
    anglplus = mirror_angle(tplus_index);
    angle(i,1) = ((anglminus-anglplus)/(tminus-tplus))*time(i) + ((tminus*anglplus)-(tplus*anglminus))/(tminus-tplus);
end

% Closing the time and position arrays
time = [mirror_time(1) time mirror_time(end)];
angle = [mirror_angle(1) ;angle ;mirror_angle(end)]; 
angle_filtered = filtfilt([1 1 1 1 1], 5,angle);

figure(1)
hold on
plot(mirror_time, mirror_angle, 'r-')
plot(time,angle,'b')
plot(time,angle_filtered,'k')
hold off
legend('raw data','interpolated data','filtered data')
xlabel('Time [ms]')
ylabel('Mirror angle [mrad]')

%% LASER SHEET POSITION

angle_filtered = angle_filtered.*1e-3;
sheet_pos_mirror = Am.*tan(angle_filtered).*sin(angle_filtered)./sin(deg2rad(45)+angle_filtered); % Initial position of the laser sheet taking account of the geometrical offset
sheet_pos_mirror = sheet_pos_mirror.*1e3; % [m] => [mm]

figure(2)
plot(time,sheet_pos_mirror);
xlabel('Time [ms]')
ylabel('Laser sheet position on the mirror [mm]')

sheet_angle = 2*angle.*1e-3; % Laser sheet angle [rad]
sheet_pos_view1 = mirr2view_dist*tan(sheet_angle); % Laser sheet position front of the viewport [m]

theta_2 = asin((n_air/n_viewport)*sin(sheet_angle));
sheet_pos_view2 = sheet_pos_view1 + view_thick*tan(theta_2); % Laser sheet position back of the viewport [m] 

theta_3 = asin((n_viewport/n_gas)*sin(theta_2));
sheet_pos_sphere = sheet_pos_view2 + sphere_rad*tan(theta_3); % Laser sheet position at the sphere center [m]

% Coefficient for the sheet local position in the sphere
pos_coeff_a = ((sheet_pos_sphere-sheet_pos_view2)*1e3)/100; % [mm/mm]
pos_coeff_b = sheet_pos_view2*1e3; % [mm]

figure(3)
legend('Interpreter','latex')
hold on
plot(time,sheet_pos_view1,'DisplayName','Front of the viewport')
plot(time,sheet_pos_view2,'DisplayName','Back of the viewport')
plot(time,sheet_pos_sphere,'DisplayName','Center of the sphere')
hold off
xlabel('Time [ms]')
ylabel('Laser sheet position [m]')

optical_axis = [0 mirr2view_dist mirr2view_dist+view_thick mirr2view_dist+view_thick+sphere_rad]*100; % Optical axis coordinates [cm] 
last_sheet = [0 sheet_pos_view1(end) sheet_pos_view2(end) sheet_pos_sphere(end)]./1e-3;
random_sheet = [0 sheet_pos_view1(end-5) sheet_pos_view2(end-5) sheet_pos_sphere(end-5)]./1e-3;

figure(4)
hold on
plot(optical_axis,last_sheet)
plot(optical_axis,random_sheet)
plot([mirr2view_dist mirr2view_dist]*100,[-20 20],'k')
plot([mirr2view_dist+view_thick mirr2view_dist+view_thick]*100,[-20 20],'k')
hold off
xlabel('Position on the setup [cm]')
ylabel('Laser sheet position [mm]')

save([WorkDir '\Positions'], 'time', 'angle', 'angle_filtered','sheet_pos_sphere','pos_coeff_a','pos_coeff_b')
disp('# POSITIONS SAVED #')