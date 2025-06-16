% Generating DEM
% Developer: Marcus Nobrega, Ph.D.
% Goal: Develop matrices to represent model inputs and B.C. for the 2D
% Boussinesq Model

%% Width Function
width_function = @(x)(50); % User define width function in terms of x

%% Input Data
L = 100;                            % Length of the hillslope [m]
angle = 5;                          % Hillslope angle in percentage
E = 0;                              % Datum reference [m]
n = 1;                              % Profile curvature parameter
omega = 0;                          % Plan Curvature Parameter
dx = 1;                             % Grid discretization [m]
dy = 1;                             % Grid discretization [m]
ksat = 1/3600;                      % Ksat [m/s]
sy = 0.3;                           % Specific yield
Aquif_Tickness = 2;                 % Aquifer tickness [m]
gw_IC = 0.2*Aquif_Tickness;         % GW Initial Depth [m]

Angle = rad2deg(atan(angle/100));   % Angle of the hillslope in deg
H = tan(pi/180*Angle)*100 + E;      % Maximum elevation [m]

%% Call DEM Function
[hillslope.DEM,hillslope.Area,hillslope.Ksat,hillslope.Sy, ...
    hillslope.GW_IC, ...
    hillslope.per,hillslope.dir_val,neu_val, ...
    hillslope.lsurf,X,Y] = hillslope_DEM(E,H,L,1,omega,dx ...
    ,dy,width_function, ...
    ksat,sy,gw_IC,Aquif_Tickness);

%% Flow Direction
% dx = X(1,2) - X(1,1);
% dy = Y(2,1) - Y(1,1);
% coord_outlet = [row_out, col_out];
% [f_dir,idx_fdir] = FlowDirection(hillslope.DEM,dx,dy,coord_outlet);
% [D_Matrix] = Find_D_Matrix(f_dir,coord_outlet,zeros(length(f_dir(:))));

%%
% [row_out,col_out] = find(hillslope.DEM == min(min(hillslope.DEM)));
% drainageMap = calculateDrainageMap(hillslope.DEM, row_out, col_out);
% surf(double(drainageMap))


x = 1:1:size(X,2);
y = min(min(Y)):1:max(max(Y));
z = hillslope.DEM;
z(isnan(z)) = nan;
z(z == -9999) = nan;
[px,py] = gradient(-z);

figure
contour(x,y,z)
hold on
quiver(x,y,px,py)
hold off
xlabel('x [m]','Interpreter','latex','FontSize',12)
ylabel('y [m]','Interpreter','latex','FontSize',12)
zlabel('z [m]','Interpreter','latex','FontSize',12)
colormap('turbo');
box on;
grid on
view(0,90);
clb = colorbar;
ylabel(clb,'Elevation [m]','Interpreter','latex','FontSize',14);

shading interp;
set(gca,'FontName','Garamond','FontSize',12)
exportgraphics(gcf,'Output/Input_Data/Gradient.png','Resolution',300)

% %% Call Correct DEM
% outlet_i = 1; outlet_j = 1;
% correct_dem = extract_upstream_cells(hillslope.DEM, outlet_i, outlet_j,1,X,Y);
% clearvars -except hillslope