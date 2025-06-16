%% Setups

%% General Setup
CMF.main('General_Setup_GW')

%% Post Processing
clear all
load 'Output/PERTH_CALIBRATION_400_particles_smaller';
% load 'Output/PERTH_CALIBRATION_no_uncertainty';
% load 'Output/30m_10m_10deg_hillslope';
Post_Processing_GW

%% Particle filter post processin
PF_Plots


