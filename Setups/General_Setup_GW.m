function [ cmf ] = General_Setup_GW()

% Read General Setup
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Inverse_Problems\Convergent_Hillslope\Drainage\Input_Data_GW_Model_convergent_drainage.xlsx';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Inverse_Problems\Convergent_Hillslope\Steady_State\Input_Data_GW_Model_convergent_steady_state';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Inverse_Problems\Convergent_Hillslope\Unsteady_State\Input_Data_GW_Model_convergent_unsteady_state';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Inverse_Problems\Flat_Hillslope\Input_Data_GW_Model_flat_drainage.xlsx';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\GitHub\2D_Boussinesq_Model\Input_Data_GW_Model_PERTH_Only_K';
dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\GitHub\2D_Boussinesq_Model\Input_Data_GW_Model_PERTH';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Input_Data_GW_Model';
% dir = 'C:\Users\marcu\OneDrive - University of Arizona\Documents\2D_Boussinesq_Model\Inverse_Problems\Uniform_Hillslope\Unsteady_State\(a)\Input_Data_GW_Model_uniform_a.xlsx';
flag_plot = 1;
[general_input,Setup_Name,BC_GW,GW_Ksat,GW_S,BC_GW_Index,GW_IC,Surf_Elevation,Recharge_BC,X,Y,SeepageFace,flag_PF] = Read_Excel_Input_GW(dir,flag_plot);
close all
%% General Settings
na = Setup_Name;                                    % run name
sO = 4;                                             % save option

% Timing Cotrol
t0 = general_input(1,1);                            % start time
tE = general_input(2,1);                            % end time
tSave = general_input(3,1):general_input(3,1):tE;                    % save times in seconds

%% Model Discretization

% Aquifer Tickness [m]
H = general_input(5,2);

%% Groundwater Model Discretization
ny = general_input(1,2);                                    % number of cells in x-direction
nx = general_input(2,2);                                    % number of cells in y-direction

dx = general_input(3,2)/nx;                         % size of cells in x-direction
dy = general_input(4,2)/ny;                        % size of cells in y-direction

%% DEM and Bedrock Elevation
lsurf = Surf_Elevation;                           % land surface

z0 = Surf_Elevation - H;    % Bedrock Elevation


%% Groundwater Model Parameters
mask = nan(ny,nx);               % active cells
mask(~isnan(GW_Ksat)) = 1;
mask = double(mask);
gwK = GW_Ksat;                   % saturated hydraulic conductivity

gwN = general_input(1,3)*ones(size(GW_Ksat));                      % V.G N parameter
gwSy = 0.3*ones(size(GW_Ksat));                                    % Now is state dependent
gwH = mask*H;                                                      % Aquifer Tickness [m]
gwPoro = general_input(2,3)*ones(size(GW_Ksat));                   % UZ Porosity [-]
gwAlpha = general_input(3,3)*ones(size(GW_Ksat));                  % Alpha of the UZ zone [1/m]
gw_k_UZ = general_input(4,3)*ones(ny,nx);                          % UZ damping factor (linear reservoir) [1/sec]
gwS_UZ = general_input(5,3)/1000*ones(ny,nx);                      % UZ storage for damping model
gw_irr_damping = general_input(6,3);                               % Damping factor in irrigation [-]

[~] = surf_plot(max(max(GW_Ksat)),1,'K','m/s',gwK,1,0,32,0.9,1,[0,90],X,Y);
exportgraphics(gcf,'Output/Input_Data/GW_Ksat.png','Resolution',600);

[~] = surf_plot(max(max(GW_S)),1,'S','-',GW_S,1,0,32,0.9,1,[0,90],X,Y);
exportgraphics(gcf,'Output/Input_Data/GW_S.png','Resolution',600);


%% Groundwater Model
Q = zeros(ny,nx,size(Recharge_BC(:,1),1));
for i = 1:size(Recharge_BC(:,1))
    Q(:,:,i) = Recharge_BC(i,2)*mask;       % m/s of recharge
end
changeTQ = Recharge_BC(:,1);                % changing times of groundwater recharge [sec]
% perimI = zeros(ny,nx);                    % cells which are boundary
% conditions

% No recharge in seepage cells (for numerical stability)
% Q(gwK == max(max(gwK)),:) = 0; % Null Q seepage cells

% Returning mask to 0
mask(isnan(mask)) = 0;

% BC GW Index
perimI = BC_GW_Index;
changeTDir = 0;                             % changing times on Dirichlet boundaries
changeTNeu = 0;                             % changing times of Neumann boundaries

perimIdir = [BC_GW(:,1)];                            % assign Dirichlet boundaries
boundDirValues = [BC_GW(:,2)];                       % Dirichlet boundary values
perimIneu = [BC_GW(:,3)];                            % assign Neumann boundaries
boundNeuValues = [BC_GW(:,4)];                       % Neumann boundary values

% create boundaries
gwboundaries = GWBoundaries(mask,Q,changeTQ,perimI,perimIdir,perimIneu,...
    changeTDir,changeTNeu,boundDirValues,boundNeuValues);

gwboundaries.SeepageFace = SeepageFace; % Adding seepage face

% wells & rivers
wells = Wells(false,0,0,0,0,0,0);
rivers = Rivers(false,0,0,0,0,0,0,0);

if max(max(GW_IC)) == 0 % No depths
    GW_IC(:,:) = 1e-6; % Add mwt depth 
end


%% Initial Conditions
IC_GW = GW_IC + z0;                             % matrix with initial GW depth + elevation [m]
% Hydrostatic Pressure for nodes in the GW table
% z0(1) is the bedrock elevation
% The total potential in the saturated zone is (z0(1) + hgw)
% The total potential in a node i, inside the saturated zone, 
% is (z0(1) + z(i)), assuming hydrostatic pressure
% Therefore, we can write
% Total Potential = z0(1) + z(i) + Head Pressure
% Head Pressure = Total Potential - z0(1) - z(i)
% Head Pressure = (z0(1) + hgw) - z0(1) - z(i)
% Head Pressure = hgw - z(i)


%% Solver and Timing
dt_GW = general_input(5,1);

if general_input(7,1) == 1
    doAdaptiveTime_GW = true;
else
    doAdaptiveTime_GW = false;
end


gwsolver = GWSolver();

%% Create Models

% Groundwater model
gwparameter = GWParameterization(gwK,gwSy,gwH,gwN,gwPoro,gwAlpha,gwS_UZ,gw_k_UZ,gw_irr_damping); % EDIT (adding Aquifer Tickness and V.G n)
gwmodel = GWModel(z0,lsurf,nx,ny,dx,dy,gwparameter,gwboundaries,wells,rivers,...
    IC_GW,gwsolver,dt_GW,doAdaptiveTime_GW,flag_PF,dir);

% coupled model
cModel = CoupledModel(gwmodel);

%% initalize framework
cmf = CMF(na,sO,tSave,t0,tE,cModel);

end

