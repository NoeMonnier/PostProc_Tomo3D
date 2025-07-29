%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  3D MESH RECONSTRUCTION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% WORKING DIRECTORY

addpath('.\Functions') % adding computing functions directory to path
WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D'); % Working directory = directory where images are stored
Original_WorkDir = WorkDir; % Keep in memory the original working directory
load([WorkDir '\Parameters']) % Load processing parameters
load([WorkDir '\Process_Parameters']) % Load parameters from contour process

%% 3D RECONSTRUCTION

for nSweep = nbSweep_start:nbSweep_end

    % WorkDir = SweepDir{nSweep}; % Change Working directory to current sweep
    % ContourDir = [WorkDir '\2.2_Contour']; % Directory where the contours data are stored
    % load([WorkDir '\Positions']); % Load the Z-position of each flame sheet
    % List_Contour = dir([ContourDir '\*.mat']); % List of all the contour files
    % Name_Contour = sortrows(char(List_Contour.name)); % Name of all the contour files
    % nbContours = size(Name_Contour,1); % Number of contours
    % contours = cell(1, nbContours); % Cell array to store XY-positions of the contour

%     for nContour = 1:nbContours
% 
%         load([ContourDir '\' Name_Contour(nContour,:)]); % Load the contour data
%         x_filt_int = real(contour_filt_int).*Magn; % Coordinates of the contour points + conversion in mm
%         y_filt_int = abs(imag(contour_filt_int)).*Magn;
%         if iscolumn(x_filt_int)
%             contours{nContour} = [x_filt_int ,y_filt_int]; % Store the data in contours cell array
%         else
%             contours{nContour} = [x_filt_int' ,y_filt_int']; 
%         end
% 
%     end
% 
%     % Conversion of contour data in a 3D Scatter plot
% 
%     XYZ = []; % All the point coordinates
%     for nContour = 1:nbContours
%         c = contours{nContour};
%         y = SweepIm_pos(nContour)*ones(size(c,1),1)*1e3;
%         XYZ = [XYZ; c(:,1), y,c(:,2)]; % Filling XYZ with all the points coordinates
%     end
% 
%     % Plotting the reconstructed shape
% 
%     shape = alphaShape(XYZ, 7); % 5 by default, to be adjusted
%     figure(nSweep)
%     plot(shape,'FaceColor','#FF7F00')
%     view(3); axis equal;
% 

    WorkDir = SweepDir{nSweep}; % Change Working directory to current sweep
    load([WorkDir '\Flame_Matrix']); % Load position and values matrices

    % Add an empty slice at the begining and the end of the flame
    empty_slice = zeros(sizey,1,sizex);
    V = cat(2,empty_slice,V);
    V = cat(2,V,empty_slice);

    closing_slice1_X = P(:,1,:,1);
    closing_slice2_X = P(:,end,:,1);
    closing_slice1_Z = P(:,1,:,3);
    closing_slice2_Z = P(:,end,:,3);
    closing_slice1_Y = -0.01+P(:,1,:,2);
    closing_slice2_Y = 0.01+P(:,end,:,2);

    % Extracting the grid coordinates and inserting the empty slices 
    X = cat(2,closing_slice1_X,P(:,:,:,1));
    X = cat(2,X,closing_slice2_X);
    Y = cat(2,closing_slice1_Y,P(:,:,:,2));
    Y = cat(2,Y,closing_slice2_Y);
    Z = cat(2,closing_slice1_Z,P(:,:,:,3));
    Z = cat(2,Z,closing_slice2_Z);


    % [X, Y, Z] = deal(P(:,:,:,1), P(:,:,:,2), P(:,:,:,3));  % Extraire les grilles
    val = 0.5;  % seuil pour l'objet
    
    fv = isosurface(X, Y, Z, V, val);
    stlwrite([WorkDir '\Unsmoothed_flame.stl'],fv);
    
    figure;
    patch(fv, 'FaceColor', '#FF7F00', 'EdgeColor', 'none');
    view(3); axis equal tight;
    camlight; lighting gouraud;
    title('Objet 3D reconstruit depuis coordonnées et valeurs');


% mask = V > 0.5;  % ou autre seuil
% 
% % Extraire les coordonnées de ces voxels
% X = P(:,:,:,1); Y = P(:,:,:,2); Z = P(:,:,:,3);
% x = X(mask); y = Y(mask); z = Z(mask);
% 
% figure;
% scatter3(x, y, z, 10, 'filled');  % couleur selon V
% view(3); axis equal tight;
% title('Nuage de points 3D depuis matrice de coordonnées');
% colorbar;
    

end



