%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CONTOUR DETECTION TOMOGRAPHY                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\Première passe de tests\NH3_Air_423K_1bar\3000 rpm\tir_1'); % Working directory = directory where images are stored
Original_WorkDir = WorkDir; % Keep in memory the original working directory
load([WorkDir '\Parameters']) % Load processing parameters

TruncDir = [WorkDir '\1_Truncation']; % Directory for truncated images
mkdir(TruncDir)
SweepsDir = [WorkDir '\2_Sweeps']; % Directory for sweep separation
mkdir(SweepsDir)
SweepDir = cell(nbSweep,1);
for i=1:nbSweep
    SweepDir{i} = [SweepsDir '\Sweep_' num2str(i,'%03.0f')]; % Directories for individual sweeps
    mkdir(SweepDir{i})
    BinDir = [SweepDir{i} '\2.1_Binary']; % Directory for binary images
    mkdir(BinDir)
    ContourDir = [SweepDir{i} '\2.2_Contour']; % Directory for contour data
    mkdir(ContourDir)
end

%% PARAMETERS

load('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\Première passe de tests\NH3_Air_423K_1bar\Positions.mat') % Load mirror and laser sheet position data

HIT_size = 40; % HIT zone size [mm]
HIT_size_pxl = floor(HIT_size/Magn); % HIT zone size [pxl]
HIT = true; % True by default, control if a flame is within the HIT zone

Trunc = true; % true = Truncate all images, true by default, put false if images have already been truncated
nbIm_start = 11; % At which image the post proc start
nbSweep_start = 7; % First Sweep to be processed

% Erosion parameters

smfilt = 7; % Size of the median filter
erosion_size = 2;
se = strel('disk',erosion_size,6); % Erosion structure

ScaleFilter = 0.4; % Scale of filtering the pixelisation noise (mm)

%% BACKGROUND IMAGE AND TRUNCATION

List_RawIm = dir([WorkDir '\*.tif']); % Gather all raw images in the directory
BACKGR = exist(BackgroundDir, 'dir'); % Checks if the background directory already exists
if BACKGR==0
    disp('#        NO BACKGROUND IMAGE FOUND        #')
    return
end

nbIm = size(List_RawIm,1); % Number of raw images
Name_RawIm = sortrows(char(List_RawIm.name)); % Name of the different images

if Trunc
    disp('### TRUNCATION ###')
    tic 
    for i = nbIm_start:nbIm
        disp(['Truncation ' Name_RawIm(i,:)]);
        RawIm = imread([WorkDir '\' Name_RawIm(i,:)]);
        BackgroundIm = imread([BackgroundDir '\' Name_RawIm(i,:)]);
        TruncIm = (imsubtract(RawIm,BackgroundIm)); % Substracting the background 
        imwrite(TruncIm,[TruncDir '\Trunc_' num2str(i,'%05.0f') '.tif']); % Saving the truncated image
    end
    disp('### TRUNCATION DONE ###')
    toc
end

%% SWEEP SEPARATION

List_TruncIm = dir([TruncDir '\*.tif']); % Gather all truncated images
Name_TruncIm = sortrows(char(List_TruncIm.name));
nbIm_Trunc = size(List_TruncIm,1); % Number of truncated images

dy = diff(position_filtered); % First derivative of the mirror position, used to find when the rotating direction change

disp('### SEPARATION ###')

%Im_counter = 1; % Count the total number of images
Sweep_counter = 1; % Count the total number of sweeps
tic
for nIm = 1:nbIm_Trunc
    disp(['Moving ' Name_TruncIm(nIm,:)]);   
    if nIm == 1
        TruncIm = imread([TruncDir '\' Name_TruncIm(nIm,:)]);
        imwrite(TruncIm,[SweepDir{Sweep_counter} '\' Name_TruncIm(nIm,:)]);
        %Im_counter = Im_counter + 1;
    else
        if dy(nIm)*dy(nIm-1)<0 % If the mirror change direction, fill an other sweep
            Sweep_counter = Sweep_counter + 1; 
        end
        TruncIm = imread([TruncDir '\' Name_TruncIm(nIm,:)]);
        imwrite(TruncIm,[SweepDir{Sweep_counter} '\' Name_TruncIm(nIm,:)]);
        %Im_counter = Im_counter + 1;
    end
end
disp('### SEPARATION DONE ###')
toc

%% MASKS

load electrodes_mask.mat % Loading default masks
load lasersheet_mask.mat

figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)])) % Show a raw image where the laser sheet is centered
h1 = drawpolygon(gca, 'Position', electrodes_pos); % Draw the mask for the electrodes
wait(h1);
disp('#   ELECTRODES MASK OK     #')
figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)])) 
h2 = drawpolygon(gca, 'Position', lasersheet_pos); % Draw the mask for the laser sheet
wait(h2); 
disp('#   LASER SHEET MASK OK    #')
figure,imshow(imread([WorkDir '\' Name_RawIm(28,:)]))
h3 = drawpoint(gca, 'Position', [sizex/2 sizey/2]); % Select the tip of the electrodes
wait(h3)
disp('#    ELECTRODES TIP OK     #')

Mask_electrodes = createMask(h1); % Create a binary mask for the electrodes
Mask_lasersheet = createMask(h2); % Create a binary mask for the laser sheet

electrodes_pos = h1.Position; % Extracting position
lasersheet_pos = h2.Position;
electrodestip_pos = h3.Position;

save('electrodes_mask.mat','electrodes_pos') % Save the created masks
save('lasersheet_mask.mat','lasersheet_pos')

close all

fig = figure('Visible','off');
imshow(imread([WorkDir '\' Name_RawIm(end,:)]))
h4 = drawcircle(gca,'Position',electrodestip_pos,'Radius',HIT_size_pxl); % Draw a circle corresponding to the HIT zone
Mask_HIT = createMask(h4); % Create a binary mask for the HIT zone

close all

%% SWEEP LOOP

for nSweep = nbSweep_start:nbSweep

    if HIT % If the flame is still inside the HIT zone

% SELECTING SUITABLE IMAGES %

        WorkDir = SweepDir{nSweep};
        
        List_SweepIm = dir([WorkDir '\*.tif']); % Gather all images in the sweep
        Name_SweepIm = sortrows(char(List_SweepIm.name));
        nbIm_Sweep = size(Name_SweepIm,1); 
        SweepIm_pos = [];
        
        for nIm=1:nbIm_Sweep
            
            figure,imshow(imread([WorkDir '\' Name_SweepIm(nIm,:)]))
            button = questdlg('Flame ?','Flame selection');
            switch button
                case 'Yes'
                    SweepIm_pos = [SweepIm_pos position(str2double(Name_SweepIm(nIm,7:end-4))-nbIm_start)];
                case 'No'
                    delete([WorkDir '\' Name_SweepIm(nIm,:)])
                case 'Cancel'
                    break
            end     
    
        end
        close all
        save([WorkDir '\Positions'],'SweepIm_pos') % Saving position of the images with flame

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

        thresholdtemp = graythresh(TruncIm_first); % Temporary threshold used for the loop
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
            BinIm_first = BinIm_first.*iminv(Mask_electrodes); 
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
            BinIm_last = BinIm_last.*iminv(Mask_electrodes); 
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
        
            figure(3)
            subplot(221); imagesc(RawIm_first); colormap gray; title('First image : raw'); axis off;
            subplot(222); imagesc(BinIm_first); title('First image : binarized'); axis off; 
            subplot(223); imagesc(RawIm_last); title('Last image : raw'); axis off; 
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
                ProcessIm = im2double(imread([WorkDir '\' Name_SweepIm(nIm,:)])); % Reading the image to process
                BinIm = imbinarize(ProcessIm,BinThreshold); % Binarize the full image with the computed threshold
                BinIm = iminv(BinIm); % Returns the negative of the image
                BinIm = BinIm.*Mask_lasersheet; % Sets pixel outside the lasersheet mask to black
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
                BinIm = BinIm.*iminv(Mask_electrodes); % Removing electrodes using the masks
    
                % Erode and Dilate the image to fill the gaps and filter the binarisation noise
            
                BinIm = imerode(BinIm,se); % Erosion and dilatation with mask se to remove binarisation noise
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
                else
                    imwrite(BinIm,strcat([WorkDir '\2.1_Binary\Binary_'],Name_SweepIm(end,7:end-4),'.TIF'),'TIF'); % If no save the image and process the next one
                end
            end
            disp('### BINARISATION DONE ###')
            toc
        end
    end

% CONTOUR DETECTION %

    

end

disp('### DONE ###')