% This code setups the Particle Filter Algorithm
% If your system is different, you might change the function
% Particle_Filter_CM() accordingly
% Developed by: Marcus Nobrega
% 10/10/2024

% Create Animation
% if flag_PF == 1 && time == 0 % Create the plots at the first time-step
%     name_plot = 'Particle_Filter_Plot';
%     obj_video = VideoWriter(['Output/Animations/' name_plot],'MPEG-4');
%     obj_video.Quality = 100;
%     obj_video.FrameRate = 10;
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
%     tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
%     open(obj_video)
% end

% Original Model Parameters
% original_K_UZ = obj.coupledModel.uzmodel.parameter.K;
% original_Alpha = obj.coupledModel.uzmodel.parameter.Alpha;
% original_N = obj.coupledModel.uzmodel.parameter.N;

% Assigning Parameters to the Model (GW)
original_K_GW = obj.coupledModel.gwmodel.parameter.K;
original_N = obj.coupledModel.gwmodel.parameter.N;
original_Poro = obj.coupledModel.gwmodel.parameter.Poro;
original_Alpha = obj.coupledModel.gwmodel.parameter.Alpha;
original_k = obj.coupledModel.gwmodel.parameter.k;

% Prior States
% original_h_UZ = obj.coupledModel.uzmodel(:).h';
% original_WC_UZ = obj.coupledModel.uzmodel(:).WC';
original_h_GW = obj.coupledModel.gwmodel.h;
original_S_UZ = obj.coupledModel.gwmodel.parameter.S_UZ;
original_Sy = obj.coupledModel.gwmodel.parameter.Sy;

if flag_PF == 1
    % Intialiazing Cluster
    % parpool(4) % See the maximum number of  workers available in your machine

    % Runs the Particle Filter of the Coupled Model
    Particle_Filter_BoussinesqModel() 

    % Plotting Frame
    f = getframe(gcf);
    % writeVideo(obj_video,f);

    % Refreshing Model Parameters to run with initial values
    obj.coupledModel.gwmodel.parameter.K = original_K_GW;
    obj.coupledModel.gwmodel.parameter.N = original_N;
    obj.coupledModel.gwmodel.parameter.Poro = original_Poro;
    obj.coupledModel.gwmodel.parameter.Alpha = original_Alpha;
    obj.coupledModel.gwmodel.parameter.k = original_k; 
    obj.coupledModel.gwmodel.parameter.Sy = original_Sy;
else
    % Parameters
    obj.coupledModel.gwmodel.parameter.K = original_K_GW;
    obj.coupledModel.gwmodel.parameter.N = original_N;
    obj.coupledModel.gwmodel.parameter.Poro = original_Poro;
    obj.coupledModel.gwmodel.parameter.Alpha = original_Alpha;
    obj.coupledModel.gwmodel.parameter.k = original_k;   
    obj.coupledModel.gwmodel.parameter.Sy = original_Sy;
end