%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     PROCESSING PARAMETERS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% PARAMETERS

WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\'); % Working directory = directory where images are stored
BackgroundDir = ('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\NH3_Air_150C_1bar_phi1\Background'); % Directory containing the background images
LaserDir = ('C:\Users\noe.monnier\Documents\Turbulent analysis\Tomo3D\NH3_Air_150C_1bar_phi1\LaserSheet'); % Directory containing laser porfile matrix

% Camera settings
F = 20000; % Camera framerate [fps]
Magn = 0.1716; % Magnification of the optical setup [mm/pxl]
sizex = 512; % Image width [pxl]
sizey = 464; % Image height [pxl]
nbSweep = 25; % Number of recorded sweeps

% Gas properties

rho_b = 0.139977008/1e6; % Density of burnt gas [g/mm^3] % NH3 ER 1 423K 1bar
rho_u = 0.746744729/1e6; % Density of unbrunt gas [g/mm^3]
sL_0 = 0.133576117; % NH3/Air 1bar 423K

fig_title = "NH$_3$/Air $\varphi$=1 P=1 bar T=423 K u'=0.33 m/s";

save([WorkDir '\Parameters'], 'BackgroundDir', 'LaserDir', 'F','Magn','sizex','sizey','nbSweep','rho_b','rho_u','sL_0','fig_title') % Save the processing parameters

disp('#  PARAMETERS SAVED  #')