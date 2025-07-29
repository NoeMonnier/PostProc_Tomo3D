%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CONTOUR DETECTION TOMOGRAPHY                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D'); % Working directory = directory where images are stored
Original_WorkDir = WorkDir; % Keep in memory the original working directory

TruncDir = [WorkDir '\1_Truncation']; % Directory for truncated images
SweepsDir = [WorkDir '\2_Sweeps']; % Directory for sweep separation

%% PARAMETERS

load([WorkDir '\Parameters']) % Load processing parameters
load([PositionsDir '\Positions.mat']) % Load laser sheet and mirror positions data

HIT_size = 40; % HIT zone size [mm]
HIT_size_pxl = floor(HIT_size/Magn); % HIT zone size [pxl]
HIT = true; % True by default, control if a flame is within the HIT zone

nbSweep_start = 10; % First Sweep to be processed

% Erosion parameters
smfilt = 7; % Size of the median filter
erosion_size = 2;
se = strel('disk',erosion_size,6); % Erosion structure

ScaleFilter = 0.4; % Scale of filtering the pixelisation noise (mm)

%% BACKGROUND IMAGE AND LASER INTENSITY PROFILE

List_RawIm = dir([WorkDir '\*.tif']); % Gather all raw images in the directory
BACKGR = exist(BackgroundDir, 'dir'); % Checks if the background directory already exists
if BACKGR==0
    disp('#        NO BACKGROUND IMAGE FOUND        #')
    return
end

nbIm = size(List_RawIm,1); % Number of raw images
Name_RawIm = sortrows(char(List_RawIm.name)); % Name of the different images

TRUNC = exist(TruncDir,'dir'); % Check if truncated image directory exists
if TRUNC==0 % If directory doesn't exist
    mkdir(TruncDir) % Create the directory
    disp('### TRUNCATION ###')
    tic 
    for i = 1:nbIm % Truncate all the image and save them in TruncDir
        disp(['Truncation ' Name_RawIm(i,:)]);
        RawIm = imread([WorkDir '\' Name_RawIm(i,:)]);
        BackgroundIm = imread([BackgroundDir '\' Name_RawIm(i,:)]);
        TruncIm = (imsubtract(RawIm,BackgroundIm)); % Substracting the background 
        imwrite(RawIm,[TruncDir '\Trunc_' num2str(i,'%05.0f') '.tif']); % Saving the truncated image
    end
    disp('### TRUNCATION DONE ###')
    toc
end

load([LaserDir '\Laser_Profile.mat'])

% LASER = exist([WorkDir '\Laser_Profile.mat'], 'file'); % Check if the laser intensity profile exists
% if LASER==0
%     disp('### LASER IMAGE FOLDER SELECTION ###')
%     LaserDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\Première passe de tests\NH3_Air_423K_1bar\3000 rpm'); % Directory where the laser images are stored
%     List_LaserIm = dir([LaserDir '\*.tif']); 
%     Name_LaserIm = sortrows(char(List_LaserIm.name)); % Name of the different images
%     nbIm_Laser = size(List_LaserIm,1); % Number of laser image
%     laserProfile = 0; % Initialising varaibles
%     laserProfile_filt = 0;
%     for nIm = 1:nbIm_Laser
%         LaserIm = im2double(imread([LaserDir '\' Name_LaserIm(nIm,:)]));
%         laserProfile = laserProfile + mean(LaserIm,2);
%         laserProfile_filt = laserProfile_filt + filtfilt(ones(1,100),100,mean(LaserIm,2));
%     end
%     laserProfile = laserProfile./10;
%     laserProfile_filt = laserProfile_filt./10;
%     save([WorkDir '\Laser_Profile.mat'],'laserProfile', 'laserProfile_filt');
% else
%     load([WorkDir '\Laser_Profile']) % Load laser profile data
% end

% Plot for visual check 
% figure(1)
% hold on
% plot(laserProfile)
% plot(laserProfile_filt,'-r')
% hold off
% legend('Laser profile','Filtered laser profile')
% saveas(gcf, [WorkDir '\Laser_Profile.fig']);

%% SWEEP SEPARATION

List_TruncIm = dir([TruncDir '\*.tif']); % Gather all truncated images
Name_TruncIm = sortrows(char(List_TruncIm.name));

dy = diff(angle_filtered); % First derivative of the mirror position, used to find when the rotating direction change

SWEEP = exist(SweepsDir,'dir'); % Check if Sweeps directory exists
if SWEEP==0
    disp('### SEPARATION ###')
    mkdir(SweepsDir) % Create diretory to store all sweeps data    
    Sweep_counter = 1; % Count the total number of sweeps
    Im_counter = 1; % Count the number of images in a sweep
    tic
    for nIm = 1:nbIm_Laser
        disp(['Moving ' Name_TruncIm(nIm,:)]);   
        if nIm == 1 % Initialising the first sweep directory and moving the first image
            SweepDir{1} = [SweepsDir '\Sweep_' num2str(Sweep_counter,'%03.0f')]; % Create the first sweep directory
            mkdir(SweepDir{1})
            CorrecDir = [SweepDir{1} '\2.1_Correc']; % Directory for laser sheet corrected images
            mkdir(CorrecDir)
            AdjustDir = [SweepDir{1} '\2.2_Adjust']; % Directory for contrast adjusted images
            mkdir(AdjustDir)
            BinDir = [SweepDir{1} '\2.3_Binary']; % Directory for binary images
            mkdir(BinDir)
            ContourDir = [SweepDir{1} '\2.4_Contour']; % Directory for contour data
            mkdir(ContourDir)
            TruncIm = imread([TruncDir '\' Name_TruncIm(nIm,:)]);
            imwrite(TruncIm,[SweepDir{Sweep_counter} '\' Name_TruncIm(nIm,:)]);
            nIm_Global(1,1) = nIm; % Store the image number in a global matrix
            Im_counter = Im_counter + 1;
        else
            if dy(nIm)*dy(nIm-1)<0 % If the mirror change direction
                Sweep_counter = Sweep_counter + 1; % Create and fill the next sweep directory
                SweepDir{Sweep_counter} = [SweepsDir '\Sweep_' num2str(Sweep_counter,'%03.0f')]; % Directories for individual sweeps
                mkdir(SweepDir{Sweep_counter})
                CorrecDir = [SweepDir{Sweep_counter} '\2.1_Correc']; % Directory for laser sheet corrected images
                mkdir(CorrecDir)
                AdjustDir = [SweepDir{Sweep_counter} '\2.2_Adjust']; % Directory for contrast adjusted images
                mkdir(AdjustDir)
                BinDir = [SweepDir{Sweep_counter} '\2.3_Binary']; % Directory for binary images
                mkdir(BinDir)
                ContourDir = [SweepDir{Sweep_counter} '\2.4_Contour']; % Directory for contour data
                mkdir(ContourDir)
                Im_counter = 1; % Reset the image counter after changing sweep directory
            end
            TruncIm = imread([TruncDir '\' Name_TruncIm(nIm,:)]);
            imwrite(TruncIm,[SweepDir{Sweep_counter} '\' Name_TruncIm(nIm,:)]);
            nIm_Global(Im_counter,Sweep_counter) = nIm; 
            Im_counter = Im_counter + 1;
        end
    end
    disp('### SEPARATION DONE ###')
    toc
    nbSweep = Sweep_counter; % Total number of sweeps
end
%% MASKS

load electrodes_mask.mat % Loading default masks
load lasersheet_mask.mat

figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)])) % Show a raw image where the laser sheet is centered
h1 = drawpolygon(gca, 'Position', electrodes_posR); % Draw the mask for the right electrodes
wait(h1);
disp('#   RIGHT ELECTRODE MASK OK   #')
h2 = drawpolygon(gca, 'Position', electrodes_posL); % Draw the mask for the left electrodes
wait(h2);
disp('#   LEFT ELECTRODE MASK OK    #')
figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)])) 
h3 = drawpolygon(gca, 'Position', lasersheet_pos); % Draw the mask for the laser sheet
wait(h3); 
disp('#   LASER SHEET MASK OK    #')
figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)]))
h4 = drawpoint(gca, 'Position', [sizex/2 sizey/2]); % Select the tip of the electrodes
wait(h4)
disp('#    ELECTRODES TIP OK     #')

Mask_electrode_R = createMask(h1); % Create a binary mask for the right electrode
Mask_electrode_L = createMask(h2); % Create a binary mask for the left electrode
Mask_lasersheet = createMask(h3); % Create a binary mask for the laser sheet

electrodes_posR = h1.Position; % Extracting position
electrodes_posL = h2.Position;
lasersheet_pos = h3.Position;
electrodestip_pos = h4.Position;

save('electrodes_mask.mat','electrodes_posR','electrodes_posL') % Save the created masks
save('lasersheet_mask.mat','lasersheet_pos')

close all

fig = figure('Visible','off');
imshow(imread([WorkDir '\' Name_RawIm(end,:)]))
h5 = drawcircle(gca,'Position',electrodestip_pos,'Radius',HIT_size_pxl); % Draw a circle corresponding to the HIT zone
Mask_HIT = createMask(h5); % Create a binary mask for the HIT zone

close all

%% SWEEP LOOP

for nSweep = nbSweep_start:nbSweep

    if HIT % If the flame is still inside the HIT zone

% SELECTING SUITABLE IMAGES %

        WorkDir = SweepDir{nSweep};
        
        List_SweepIm = dir([WorkDir '\*.tif']); % Gather all images in the sweep
        Name_SweepIm = sortrows(char(List_SweepIm.name));
        nbIm_Sweep = size(Name_SweepIm,1); 
        SweepIm_coeff = [];
        
        for nIm=1:nbIm_Sweep
            
            figure,imshow(imadjust(imread([WorkDir '\' Name_SweepIm(nIm,:)]))) % Manually check if the image contains a flame
            button = questdlg('Flame ?','Flame selection');
            switch button
                case 'Yes'
                    SweepIm_coeff = [SweepIm_coeff; pos_coeff_a(str2double(Name_SweepIm(nIm,7:end-4))) pos_coeff_b(str2double(Name_SweepIm(nIm,7:end-4)))]; % Store the physicak position of the laser sheet
                case 'No'
                    delete([WorkDir '\' Name_SweepIm(nIm,:)])
                case 'Cancel'
                    break
            end     
    
        end
        close all
        save([WorkDir '\Positions'],'SweepIm_coeff') % Saving position of the images with flame

        % Re evaluate the image list in the directory
        List_SweepIm = dir([WorkDir '\*.tif']); % Gather all images in the sweep
        Name_SweepIm = sortrows(char(List_SweepIm.name));
        nbIm_Sweep = size(Name_SweepIm,1);

% BINARISATION THRESHOLD % 

        threshold_test = true; % Controls the loop
        coeff = 1.0; % Coefficient to adjust the threshold

        RawIm_first = imread([Original_WorkDir '\' Name_RawIm(str2double(Name_SweepIm(1,7:end-4)),:)]);
        RawIm_last = imread([Original_WorkDir '\' Name_RawIm(str2double(Name_SweepIm(end,7:end-4)),:)]);
        SweepIm_first = im2double(imread([WorkDir '\' Name_SweepIm(1,:)]));
        SweepIm_last = im2double(imread([WorkDir '\' Name_SweepIm(end,:)]));

        % Normalisation to take account of the non homogeneity of the laser sheet
        for nColl = 1:sizey
                SweepIm_first(:,nColl) = SweepIm_first(:,nColl)./(laserProfile_filt(:,nIm_Global(1,nSweep))/max(laserProfile_filt(:,nIm_Global(1,nSweep)))); 
                SweepIm_last(:,nColl) = SweepIm_last(:,nColl)./(laserProfile_filt(:,nIm_Global(1,nSweep))/max(laserProfile_filt(:,nIm_Global(1,nSweep)))); 
        end

        SweepIm_first = imadjust(SweepIm_first);
        SweepIm_last = imadjust(SweepIm_last);
        thresholdtemp = graythresh(SweepIm_first(find(Mask_lasersheet==1))); % Temporary threshold used for the loop
        fprintf('Binary threshold = %0.03f', thresholdtemp); % Display the value of the temporary threshold

        while threshold_test

            thresholdtemp = thresholdtemp*coeff; % Updating the threshold
        
            % Full binarisation process for the first and last images
            BinIm_first = imbinarize(SweepIm_first, thresholdtemp);
            BinIm_first = iminv(BinIm_first); 
            BinIm_first = BinIm_first.*Mask_lasersheet;
            % if continuous_flame
            BinIm_first = Fct_Struct_Max(BinIm_first); 
            % else
            %     BinIm_first_label = bwlabel(BinIm_first,4);
            %     nb_struct = max(BinIm_first_label(:));
            %     for num_struct = 1:nb_struct
            %         Pix_in_struct = find(BinIm_first_label==num_struct);
            %         if length(Pix_in_struct)<200
            %             BinIm_first(Pix_in_struct)=0;
            %         end
            %     end
            % end
            BinIm_first = BinIm_first.*iminv(Mask_electrode_R); 
            BinIm_first = BinIm_first.*iminv(Mask_electrode_L);
            BinIm_first = imerode(BinIm_first,se); 
            BinIm_first = imdilate(BinIm_first,se);
            %BinIm_first = Fct_Struct_Max(BinIm_first); 
            BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]);
            BinIm_first = imdilate(BinIm_first,se);
            BinIm_first = imdilate(BinIm_first,se);
            BinIm_first = imdilate(BinIm_first,se);
            BinIm_first = imdilate(BinIm_first,se);
            BinIm_first = imdilate(BinIm_first,se);
            BinIm_first = imerode(BinIm_first,se);
            BinIm_first = imerode(BinIm_first,se);
            BinIm_first = imerode(BinIm_first,se);
            BinIm_first = imerode(BinIm_first,se);
            BinIm_first = imerode(BinIm_first,se);
            BinIm_first = imfill(BinIm_first,'holes'); 
            BinIm_first = Fct_Struct_Max(BinIm_first); 
            BinIm_first = medfilt2(BinIm_first, [smfilt smfilt]); 
        
            BinIm_last = imbinarize(SweepIm_last, thresholdtemp);
            BinIm_last = iminv(BinIm_last); 
            BinIm_last = BinIm_last.*Mask_lasersheet;
            % if continuous_flame
            BinIm_last = Fct_Struct_Max(BinIm_last); 
            % else
            %     BinIm_last_label = bwlabel(BinIm_last,4);
            %     nb_struct = max(BinIm_last_label(:));
            %     for num_struct = 1:nb_struct
            %         Pix_in_struct = find(BinIm_last_label==num_struct);
            %         if length(Pix_in_struct)<200
            %             BinIm_last(Pix_in_struct)=0;
            %         end
            %     end
            % end
            BinIm_last = BinIm_last.*iminv(Mask_electrode_R);
            BinIm_last = BinIm_last.*iminv(Mask_electrode_L);
            BinIm_last = imerode(BinIm_last,se); 
            BinIm_last = imdilate(BinIm_last,se);
            %BinIm_last = Fct_Struct_Max(BinIm_last); 
            BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]);
            BinIm_last = imdilate(BinIm_last,se);
            BinIm_last = imdilate(BinIm_last,se);
            BinIm_last = imdilate(BinIm_last,se);
            BinIm_last = imdilate(BinIm_last,se);
            BinIm_last = imdilate(BinIm_last,se);
            BinIm_last = imerode(BinIm_last,se);
            BinIm_last = imerode(BinIm_last,se);
            BinIm_last = imerode(BinIm_last,se);
            BinIm_last = imerode(BinIm_last,se);
            BinIm_last = imerode(BinIm_last,se);
            BinIm_last = imfill(BinIm_last,'holes'); 
            BinIm_last = Fct_Struct_Max(BinIm_last); 
            BinIm_last = medfilt2(BinIm_last, [smfilt smfilt]); 
            
            % Contour detection for the first and last image to check if binarisation threshold is good
            fig_cont_first = figure('Visible','off');
            cont_first = imcontour(BinIm_first,1); % Contour detection
            close(fig_cont_first)

            fig_cont_last = figure('Visible','off');
            cont_last = imcontour(BinIm_last,1); % Contour detection
            close(fig_cont_last)

            x_cont_first = cont_first(1,2:end); % Contour coordinates
            y_cont_first = cont_first(2,2:end);

            x_cont_last = cont_last(1,2:end); % Contour coordinates
            y_cont_last = cont_last(2,2:end);

            % Check if contour is closed
            if x_cont_first(end)==x_cont_first(1) && y_cont_first(end) == y_cont_first(1)
                if x_cont_first(1)== x_cont_first(end) % If contour extremities are aligned
                    x_cont_first = [x_cont_first x_cont_first(1)]; % closing manually the contour
                    y_cont_first = [y_cont_first y_cont_first(1)];
                else
                    slope = (y_cont_first(end)-y_cont_first(1))/(x_cont_first(end)-x_cont_first(1)); % Closing the contour using linear segment
                    y_intersect = y_cont_first(1)-slope*x_cont_first(1);
                    x_line = x_cont_first(1)+0.5:0.5:x_cont_first(end)-0.5;
                    y_line = slope.*x_line + y_intersect;
                    x_cont_first = [x_cont_first fliplr(x_line) x_cont_first(1)];
                    y_cont_first = [y_cont_first fliplr(y_line) y_cont_first(1)];
                end
            end

            if x_cont_last(end)==x_cont_last(1) && y_cont_last(end) == y_cont_last(1)
                if x_cont_last(1)== x_cont_last(end) % If contour extremities are aligned
                    x_cont_last = [x_cont_last x_cont_last(1)]; % closing manually the contour
                    y_cont_last = [y_cont_last y_cont_last(1)];
                else
                    slope = (y_cont_last(end)-y_cont_last(1))/(x_cont_last(end)-x_cont_last(1)); % Closing the contour using linear segment
                    y_intersect = y_cont_last(1)-slope*x_cont_last(1);
                    x_line = x_cont_last(1)+0.5:0.5:x_cont_last(end)-0.5;
                    y_line = slope.*x_line + y_intersect;
                    x_cont_last = [x_cont_last fliplr(x_line) x_cont_last(1)];
                    y_cont_last = [y_cont_last fliplr(y_line) y_cont_last(1)];
                end
            end

            % Remove points in the electrodes mask from the contour
            for contour_point = 1:size(x_cont_first,2)
                % Right electrode
                [in, on] = inpolygon(x_cont_first(contour_point),y_cont_first(contour_point),electrodes_posR(:,1),abs(electrodes_posR(:,2)));
                if in || on 
                    x_cont_first(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                    y_cont_first(contour_point) = NaN;
                end
                % Left electrode
                [in, on] = inpolygon(x_cont_first(contour_point),y_cont_first(contour_point),electrodes_posL(:,1),abs(electrodes_posL(:,2)));
                if in || on 
                    x_cont_first(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                    y_cont_first(contour_point) = NaN;
                end
            end

            % Extract all the non-NaN values from the contour coordinates
            x_cont_first = x_cont_first(1,~isnan(x_cont_first(1,:)));
            y_cont_first = y_cont_first(1,~isnan(y_cont_first(1,:)));

            for contour_point = 1:size(x_cont_last,2)
                % Right electrode
                [in, on] = inpolygon(x_cont_last(contour_point),y_cont_last(contour_point),electrodes_posR(:,1),abs(electrodes_posR(:,2)));
                if in || on 
                    x_cont_last(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                    y_cont_last(contour_point) = NaN;
                end
                % Left electrode
                [in, on] = inpolygon(x_cont_last(contour_point),y_cont_last(contour_point),electrodes_posL(:,1),abs(electrodes_posL(:,2)));
                if in || on 
                    x_cont_last(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                    y_cont_last(contour_point) = NaN;
                end
            end

            % Extract all the non-NaN values from the contour coordinates
            x_cont_last = x_cont_last(1,~isnan(x_cont_last(1,:)));
            y_cont_last = y_cont_last(1,~isnan(y_cont_last(1,:)));

            contour_first = x_cont_first -1i * y_cont_first; % Reconstructing the contour geometry
            contour_last = x_cont_last -1i * y_cont_last;

            % Contour filtering
            contour_filt_first = Fct_Contour_Filter(contour_first, ScaleFilter, Magn);
            contour_filt_last = Fct_Contour_Filter(contour_last, ScaleFilter, Magn);
        
            % Contour interpolation
            % Compute the distance between two successive points 
            % If interpolation is good the distance should be around 0.5 
            % Else interpolate again
            contour_filt_int_first = contour_filt_first;
            contour_filt_int_last = contour_filt_last;

            % while min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.9
            while ((max(sqrt(diff(real(contour_filt_int_first)).^2+diff(imag(contour_filt_int_first)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int_first)).^2+diff(imag(contour_filt_int_first)).^2))<0.45))
                contour_temp = contour_filt_int_first(:)';
                contour_filt_int_first = Fct_ContourInterpSpline(contour_temp,1,1);
            end

            while ((max(sqrt(diff(real(contour_filt_int_last)).^2+diff(imag(contour_filt_int_last)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int_last)).^2+diff(imag(contour_filt_int_last)).^2))<0.45))
                contour_temp = contour_filt_int_last(:)';
                contour_filt_int_last = Fct_ContourInterpSpline(contour_temp,1,1);
            end
        
            % Extract the filtered contour points coordinates
            x_filt_int_first = real(contour_filt_int_first);
            y_filt_int_first = abs(imag(contour_filt_int_first));
            x_filt_int_last = real(contour_filt_int_last);
            y_filt_int_last = abs(imag(contour_filt_int_last));
        
            figure(3)
            subplot(221); imagesc(imadjust(RawIm_first)); colormap gray; title('First image : raw'); hold on; plot(x_filt_int_first,y_filt_int_first,'r-'); hold off; axis off
            subplot(222); imagesc(BinIm_first); title('First image : binarized'); axis off; 
            subplot(223); imagesc(imadjust(RawIm_last)); colormap gray; title('Last image : raw'); hold on; plot(x_filt_int_last,y_filt_int_last,'r-'); hold off; axis off
            subplot(224); imagesc(BinIm_last); title('Last image : binarized'); axis off; 
            pause(3)
        
            button = questdlg('Threshold OK ?','BinThreshold');
                switch button
                    case 'Yes'
                        BinThreshold = thresholdtemp;
                        threshold_test = false;
                    case 'No'
                        Para = inputdlg({'New coefficient : '}, 'coeff', 1.0, {'1.0'});
                        coeff = str2double(Para{1});
                    case 'Cancel'
                        break
                end     
        end
        close(3)

% BINARISATION PROCESS %

        disp('# BINARISATION #')

        tic
        for nIm = 1:nbIm_Sweep
            if HIT 
                disp(['Binarisation ' Name_SweepIm(nIm,:)]);
                BinIm = im2double(imread([WorkDir '\' Name_SweepIm(nIm,:)])); % Reading the image to process
                % Normalisation to take account of the non homogeneity of the laser sheet
                for nColl = 1:sizey
                    BinIm(:,nColl) = BinIm(:,nColl)./(laserProfile_filt(:,nIm_Global(1,nSweep))/max(laserProfile_filt(:,nIm_Global(1,nSweep)))); 
                end
                imwrite(BinIm,strcat([WorkDir '\2.1_Correc\Correc_'],Name_SweepIm(nIm,7:end-4),'.tif'),'TIF'); % save the image
                BinIm = imadjust(BinIm); % Increase the contrast of the image for easier binarization
                imwrite(BinIm,strcat([WorkDir '\2.2_Adjust\Adjust_'],Name_SweepIm(nIm,7:end-4),'.tif'),'TIF'); % If no save the image and process the next one
                BinIm = imbinarize(BinIm,BinThreshold); % Binarize the full image with the computed threshold
                BinIm = iminv(BinIm); % Returns the negative of the image
                % BinIm = BinIm.*Mask_lasersheet; % Sets pixel outside the lasersheet mask to black
                % if continuous_flame
                    BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure in the image
                % else
                %     BinIm_label = bwlabel(BinIm,4);
                %     nb_struct = max(BinIm_label(:));
                %     for num_struct = 1:nb_struct
                %         Pix_in_struct = find(BinIm_label==num_struct);
                %         if length(Pix_in_struct)<200
                %             BinIm(Pix_in_struct)=0;
                %         end
                %     end
                % end
                BinIm = BinIm.*iminv(Mask_electrode_R); % Removing electrodes using the masks
                BinIm = BinIm.*iminv(Mask_electrode_L);
    
                % Erode and Dilate the image to fill the gaps and filter the binarisation noise
            
                BinIm = imerode(BinIm,se); % Erosion and dilatation with structure se to remove binarisation noise
                BinIm = imdilate(BinIm,se);
                % BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure in the image
                BinIm = medfilt2(BinIm, [smfilt smfilt]); % Median filter to remove the last binarisation noise
            
                BinIm = imdilate(BinIm,se);
                BinIm = imdilate(BinIm,se);
                BinIm = imdilate(BinIm,se);
                BinIm = imdilate(BinIm,se);
                BinIm = imdilate(BinIm,se);
                BinIm = imerode(BinIm,se);
                BinIm = imerode(BinIm,se);
                BinIm = imerode(BinIm,se);
                BinIm = imerode(BinIm,se);
                BinIm = imerode(BinIm,se);
            
                BinIm = imfill(BinIm,'holes'); % Fill the detection gaps in the image
                BinIm = medfilt2(BinIm, [5 5]); % Median filter to remove the last binarisation noise
                BinIm = Fct_Struct_Max(BinIm); % Keeping the largest structure on the image (normally the flame)
        
                % Check if the detected flame exit the HIT zone
                if max(BinIm.*(1-Mask_HIT),[],'all')==1
                    HIT = false; % If yes don't save and don't process aditionnal images
                    disp('Flame exited the HIT zone')
                    nbSweep_end = nSweep;
                else
                    imwrite(BinIm,strcat([WorkDir '\2.3_Binary\Binary_'],Name_SweepIm(nIm,7:end-4),'.tif'),'TIF'); % If no save the image and process the next one
                end
            end
        end
        disp('### BINARISATION DONE ###')
        toc
    end

% CONTOUR DETECTION %

    disp("### CONTOUR DETECTION ###")
    
    WorkDir = [SweepDir{nSweep} '\2.3_Binary'];
    ContourDir = [SweepDir{nSweep} '\2.4_Contour'];
    List_BinIm = dir([WorkDir '\*.tif']); % Gather all binary images from the directory
    Name_BinIm = sortrows(char(List_BinIm.name)); % Name of all binary images
    nbIm_Bin = size(List_BinIm,1); % Number of images to process

    tic
    for nIm=1:nbIm_Bin
        disp(['Contour detection ' Name_BinIm(nIm,:)]);
        BinIm = im2double(imread([WorkDir '\' Name_BinIm(nIm,:)])); % Reading the binary image
        
        fig_cont = figure('Visible','off');
        cont = imcontour(BinIm,1); % Contour detection
        close(fig_cont)
    
        x_cont = cont(1,2:end); % Contour coordinates
        y_cont = cont(2,2:end);

         % Check if contour is closed
        if x_cont(end)==x_cont(1) && y_cont(end) == y_cont(1)
            if x_cont(1)== x_cont(end) % If contour extremities are aligned
                x_cont = [x_cont x_cont(1)]; % closing manually the contour
                y_cont = [y_cont y_cont(1)];
            else
                slope = (y_cont(end)-y_cont(1))/(x_cont(end)-x_cont(1)); % Closing the contour using linear segment
                y_intersect = y_cont(1)-slope*x_cont(1);
                x_line = x_cont(1)+0.5:0.5:x_cont(end)-0.5;
                y_line = slope.*x_line + y_intersect;
                x_cont = [x_cont fliplr(x_line) x_cont(1)];
                y_cont = [y_cont fliplr(y_line) y_cont(1)];
            end
        end

        % Remove points in the electrodes mask from the contour
        for contour_point = 1:size(x_cont,2)
            % Right electrode
            [in, on] = inpolygon(x_cont(contour_point),y_cont(contour_point),electrodes_posR(:,1),abs(electrodes_posR(:,2)));
            if in || on 
                x_cont(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                y_cont(contour_point) = NaN;
            end
            % Left electrode
            [in, on] = inpolygon(x_cont(contour_point),y_cont(contour_point),electrodes_posL(:,1),abs(electrodes_posL(:,2)));
            if in || on 
                x_cont(contour_point) = NaN; % replacing values by NaN to keep the array size consistent
                y_cont(contour_point) = NaN;
            end
        end

        % Extract all the non-NaN values from the contour coordinates
        x_cont = x_cont(1,~isnan(x_cont(1,:)));
        y_cont = y_cont(1,~isnan(y_cont(1,:)));

        contour = x_cont -1i * y_cont; % Reconstructing the contour geometry

        % Contour filtering
        contour_filt = Fct_Contour_Filter(contour, ScaleFilter, Magn);
    
        % Contour interpolation
        % Compute the distance between two successive points 
        % If interpolation is good the distance should be around 0.5 
        % Else interpolate again
        contour_filt_int = contour_filt;
        % while min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.9
        while ((max(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))>0.55) && (min(sqrt(diff(real(contour_filt_int)).^2+diff(imag(contour_filt_int)).^2))<0.45))
            contour_temp = contour_filt_int(:)';
            contour_filt_int = Fct_ContourInterpSpline(contour_temp,1,1);
        end
    
        % Extract the filtered contour points coordinates
        x_filt_int = real(contour_filt_int);
        y_filt_int = abs(imag(contour_filt_int));

        % Plot the contour on the corresponding raw image for validation
        RawIm = imread([Original_WorkDir '\tir' Name_BinIm(nIm,8:end-4) '.TIF']);
        contour_plot = figure('Visible','off');
        imagesc(RawIm); axis equal; title(['Image n°' Name_BinIm(nIm,7:end-4)]); colormap gray; hold on; plot(x_filt_int,y_filt_int,'r-');
        saveas(contour_plot,[SweepDir{nSweep} '\2.4_Contour\Image' Name_BinIm(nIm,7:end-4) '.tif']);
        close(contour_plot)
    
        save([ContourDir '\Contour_' Name_BinIm(nIm,7:end-4)], 'contour','contour_filt','contour_filt_int');

    end
    disp('### CONTOUR DETECTION DONE ###')
    toc

% V and P Matrix

    V = zeros(sizey,nbIm_Bin,sizex); % Value matrix : 1 if flame 0 if no flame
    P = zeros(sizey,nbIm_Bin,sizex,3); % Position matrix : store the XYZ position of each V points
    
    for nIm=1:nbIm_Bin
        BinIm = imread([WorkDir '\' Name_BinIm(nIm,:)]); % Reading the binary image
        for i = 1:sizey
            for j = 1:sizex
                V(i,nIm,j) = BinIm(i,j); % Store the slice in V matrix
            end
        end
        x = ((1:sizex)-round(electrodestip_pos(1)))*Magn; % X coordinates of each point in the image [mm]
        z = ((1:sizey)-round(electrodestip_pos(2)))*Magn; % Z coordinates
        % Store the coordinates in P
        for i=1:sizey
            P(i,nIm,:,1) = x;
        end

        for i = 1:sizex
            P(:,nIm,i,2) = SweepIm_coeff(nIm,1).*P(:,nIm,i,1) + SweepIm_coeff(nIm,2); % Compute the Y position during propagation through the sphere
        end

        for i=1:sizex
            P(:,nIm,i,3) = z;
        end
        

    end

    save([SweepDir{nSweep} '\Flame_Matrix'], 'V', 'P')

end

save([Original_WorkDir '\Process_parameters'], 'nbSweep', 'nbSweep_start', 'nbSweep_end', 'SweepDir');

disp('### DONE ###')

