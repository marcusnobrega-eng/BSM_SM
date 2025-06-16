%% GW Model

% 1.0 - Read Excel Input Data
dir = 'C:\Users\Marcus\Desktop\Desktop_Folder\LEO\h3d-brandhorst-erdal-main-CoupledIterative\CoupledIterative\GW_Model_Only\Input_Data_GW_Model.xlsx';
[general_input,Setup_Name,BC_GW,GW_Ksat,GW_S,BC_GW_Index,GW_IC,Surf_Elevation,Irrigation_BC] = Read_Excel_Input_GW(dir);

% 2.0 - Process