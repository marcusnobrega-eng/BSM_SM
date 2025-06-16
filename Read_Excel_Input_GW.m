function [general_input,Setup_Name,BC_GW,GW_Ksat,GW_S,BC_GW_Index,GW_IC,Surf_Elevation,Recharge_BC,X,Y,SeepageFace,flag_PF] = Read_Excel_Input_GW(dir,flag_plot)

% Read Input Data from Directory
% dir = 'C:\Users\Marcus\Desktop\Desktop_Folder\LEO\h3d-brandhorst-erdal-main-CoupledIterative\CoupledIterative\Input_Data_Coupled_Model.xlsx';
input_data = readtable(dir,'Sheet','Input_Data'); % Reading Input data from Excel
BC_data = readtable(dir,'Sheet','Boundary_Conditions'); % Reading B.C
GW_depth = readtable(dir,'Sheet','GW_Depth');
GW_Ksat = readtable(dir,'Sheet','GW_Ksat');
GW_S = readtable(dir,'Sheet','GW_S');
DEM_data = readtable(dir,'Sheet','DEM');
Recharge_data = readtable(dir,'Sheet','Inflow_BC');

try 
    mkdir Output/Input_Data
catch 
    warning('Folder (Input_Data) already created')
end

%% General Input
% Organize States of each Input Data
general_input = zeros(7,4);
% Col 1 - Time-Stepping
general_input(1:7,1) = table2array(input_data(2:8,2));
% % Col 2 - GW Discretization and Parameters
general_input(1:5,2) = table2array(input_data(2:6,5));
% Col 3 - UZ Surrogate Model
general_input(1:6,3) = table2array(input_data(8:13,5));
% Setup Name
Setup_Name = char(table2array(input_data(2,8)));
% flag_PF
flag_PF = table2array(input_data(7,2));

%% Boundary Conditions
% B.C Indexes and Values
temp = table2array(BC_data(:,1:4));
n_not_nan = max(sum(~isnan(temp(:,1))), sum(~isnan(temp(:,3))));
temp((n_not_nan+1):end,:) = [];
% Dirichlet Index
n_dirc = length(temp(:,1)) - sum(isnan(temp(:,1)));
BC_GW(1:n_dirc,1) = temp(1:n_dirc,1);
% Dirchlet Values
BC_GW(1:n_dirc,2) = temp(1:n_dirc,2);
% Neumann Index
n_neu = length(temp(:,3)) - sum(isnan(temp(:,3)));
BC_GW(1:n_neu,3) = temp(1:n_neu,3);
BC_GW(1:n_neu,4) = temp(1:n_neu,4);
% Matrix Values
n_rows = size(BC_data,1);
n_cols = size(BC_data,2) - 6;
BC_GW_Index = table2array(BC_data(1:n_rows,(7:(7+n_cols-1))));
SeepageFace = zeros(size(BC_GW_Index));
for i = 1:n_dirc
    [row,col] = find(BC_GW_Index == BC_GW(i,1));
    SeepageFace(row,col) = 1;
end
clear temp

%% Groundwater Initial Depth
temp = table2array(GW_depth(:,2:end));
GW_IC = temp;
clear temp

%% Groundwater Ksat
temp = table2array(GW_Ksat(:,2:end));
GW_Ksat = temp;
clear temp

%% Groundwater Specific Yield
temp = table2array(GW_S(:,2:end));
GW_S = temp;
clear temp
%% DEM
temp = table2array(DEM_data(:,2:end));
Surf_Elevation = temp;

%% Treating Nan Values
idx_nan = Surf_Elevation == -9999;
Surf_Elevation(idx_nan) = nan;
GW_S(idx_nan) = nan;
GW_Ksat(idx_nan) = nan;
BC_GW_Index(idx_nan) = nan;
GW_IC(idx_nan) = nan;

[X,Y] = meshgrid(1:1:size(BC_GW_Index,2),1:1:size(BC_GW_Index,1));
X = X*(general_input(3,2)/(size(X,2))); Y = Y*(general_input(4,2)/(size(X,1)));
if flag_plot == 1
surf_plot(max(max(Surf_Elevation)),1,'\mathrm{DEM}','m',Surf_Elevation,0,0,256,1,1,[0,90],X,Y);
exportgraphics(gcf,'Output/Input_Data/DEM.png','Resolution',600);
end

%% Inflow B.C
temp = table2array(Recharge_data);
net_rainfall = -(temp(:,2) - temp(:,3))/1000/3600; % mm/h to m/s
Recharge_BC(:,1) = temp(:,1)*60; % Minutes to Seconds
Recharge_BC(:,2) = net_rainfall;

end