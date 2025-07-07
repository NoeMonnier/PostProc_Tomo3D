close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('./Functions/') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis'); % Working directory, where the laser swing videos are located
Original_WorkDir = WorkDir; % Keep in memory the original working directory

%% PARAMETERS

BackgroundDir = ('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\NH3_Air_150C_1bar_phi1\Background'); % Directory containing the background images
nb_shots = 10; % Number of experimental shots
sizey = 464; % Images height

%% MEMORY ALLOCATION

nb_im=NaN(1,nb_shots); % Number of processed images for each shot

%% TRUNCATION
% All the images are truncated with a background image to remove the laser
% reflections on the elecrodes and vessel walls

BACKGR = exist(BackgroundDir, 'dir'); % Checks if the background directory already exists
if BACKGR==0
    disp('#        NO BACKGROUND IMAGE FOUND        #')
    return
end

for n_shot = 1:nb_shots

    disp(['Truncation Tir ' num2str(n_shot)]);

    LaserDir = [WorkDir '\tir_' num2str(n_shot)]; % Directory of the n-th experimental shot
    TruncDir = [LaserDir '\1-Truncation']; % Directory to store truncated images

    List_RawIm = dir([LaserDir '\*.tif']); % Gather all the images of the n-th shot
    Name_RawIm = sortrows(char(List_RawIm.name)); % Name of the gathered images
    nb_im(n_shot) = size(List_RawIm,1); % Number of raw images  

    TRUNC = exist(TruncDir, 'dir'); % Check if images have already been truncated
    if TRUNC == 0 % If images haven't already been truncated create the folder and truncate all the images

        mkdir(TruncDir)
        
        for n_im = 1:nb_im(n_shot)     

            RawIm = imread([LaserDir '\' Name_RawIm(n_im,:)]);
            BackgroundIm = imread([BackgroundDir '\' Name_RawIm(n_im,:)]);
            TruncIm = (imsubtract(RawIm,BackgroundIm)); % Substracting the background 
            imwrite(TruncIm,[TruncDir '\Trunc_' num2str(n_im,'%05.0f') '.tif']); % Saving the truncated image 

        end
    else
        disp('Images already truncated')
    end
end


%% MEAN PROFILE DETERMINATION 

nbIm_Laser = min(nb_im); %  
laserProfile = zeros(sizey,nbIm_Laser); 
laserProfile_filt = zeros(sizey,nbIm_Laser);

for n_im = 1:nbIm_Laser

    disp(['Laser profile image : ' num2str(n_im)])
    
    for n_shot = 1:nb_shots
        
        TruncDir = [WorkDir '\tir_' num2str(n_shot) '\1-Truncation']; % Directory with the truncated images
        List_LaserIm = dir([TruncDir '\*.tif']); % Gather all the images
        Name_LaserIm = sortrows(char(List_LaserIm.name)); % Name of the gathered images
        LaserIm = im2double(imread([TruncDir '\' Name_LaserIm(n_im,:)]));
        laserProfile(:,n_im) = laserProfile(:,n_im) + mean(LaserIm,2);
        laserProfile_filt(:,n_im) = laserProfile_filt(:,n_im) + filtfilt(ones(1,100),100,mean(LaserIm,2));
    end

end

laserProfile = laserProfile./10;
laserProfile_filt = laserProfile_filt./10;
save([WorkDir '\Laser_Profile.mat'],'nbIm_Laser','laserProfile', 'laserProfile_filt');

disp('DONE')