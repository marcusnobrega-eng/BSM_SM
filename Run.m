%% Setups

%% General Setup
CMF.main('General_Setup_GW')

%% Post Processing
clear all
load 'Output/PERTH_CALIBRATION_400_particles_smaller';
Post_Processing_GW

%% Particle filter post processinG
PF_Plots


