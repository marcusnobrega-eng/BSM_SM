% Post-Processing Code
% Developer: Marcus Nobrega
%% 1.0 Load Data
clear all
% load 'Output'/LEO_test.mat
% load 'Output'/Bucket.mat
% load 'Output'/miniLEO.mat
% load 'Output'/LEO.mat
% load 'Output'/LEO_steady_state_05m.mat
load 'Output'/LEO_steady_state_irr_2mmhr_1day.mat

close all

%% Coloramp
[Spectrum,depth_ramp,terrain_ramp] = coloramps(); % Run coloramp function

%% 2.0 Data
time_step = tSave(1,1);                 % Seconds
time_plot_days = tSave/(3600*24);       % days

%% 3.0 Reshaping Data to 2D matrices
nGW_states = size(pressureHead,1) - size(waterContent,1);                        % Number of GW states
nUZ_states = size(pressureHead,1) - nGW_states;
ny = size(z0,1);
nx = size(z0,2);
[X,Y] = meshgrid(1:1:nx,1:1:ny);
dz_cells = cumsum(deltaCS);
nz = length(dz_cells);
zonation = length(xi);

% Pressure Head
head_pressure = pressureHead(1:nUZ_states,:);
head_pressure = reshape(head_pressure,[nz,length(unique(yi)),length(unique(xi)),size(head_pressure,2)]); % nUZ x nx x ny x T fields of pressure

% Potential Head
potential_head = pressureHead(1:nUZ_states,:);
for i = 1:length(xi)
    z_zone(i,1) = z0(yi(i),xi(i));
end
potential_head = reshape(potential_head,[nz,length(unique(yi)),length(unique(xi)),size(potential_head,2)]); % nUZ x nx x ny x T fields of pressure
z_Zone(1,:,:,:) = repmat(reshape(z_zone,[size(potential_head,2),size(potential_head,3)]),[1,1,size(potential_head,4)]);
potential_head = potential_head + z_Zone;
% Moisture
waterContent = waterContent(1:nUZ_states,:);
waterContent = reshape(waterContent,[nz,length(unique(yi)),length(unique(xi)),size(waterContent,2)]); % nz x ny x nx x T fields of pressure

% Groundwater Levels
GW_table = GW + z0; % m
GW_table_original = GW_table;
GW_Storage = squeeze(sum(sum(GW)))*dx*dy; % m3
GW_table = reshape(GW_table,[1,size(GW_table,1),size(GW_table,2),size(GW_table,3)]);
GW_original = GW;
GW = reshape(GW,[1,size(GW,1),size(GW,2),size(GW,3)]);

%% 4.9 Profile Plot
xi_unique = unique(xi);
yi_unique = unique(yi);
if size(xi_unique,1) > 1
    zone_coordinates = [yi, xi];
else
    zone_coordinates = [yi', xi'];
end

% Finding Closest Cells to the Middle of the Domain
% x
mid_x = round(nx/2);

id_plot = [57,14]; % Cells
cell_plot = [find(yi_unique == id_plot(1),1,'first') ; find(xi_unique == id_plot(2),1,'first')]';
if length(cell_plot) < 2 || isempty(cell_plot)
    warning('Choose another cell and make sure it has a UZ model associated with.')
    cell_plot = [length(yi_unique);length(xi_unique)];
end
cell_id = [yi_unique(cell_plot(1,1)), xi_unique(cell_plot(1,2))];
%% 3.0 Plots

% (a) GW

% (b) Sy

% (c) Pressure

%% 4.0 Colors
n_plots = 10;
colors = linspecer(n_plots,1);

%% Bedrock Elevation
elevation_surf(time_plot_days,'Bedrock Elevation','z','m',z0)

%% Piezometric Head Animation
plot_animation(deltaCS,H,nz,cell_plot,time_plot_days,'Potential_Data','\mathrm{Piezometric~Head}','m',potential_head,colors(1,:),cell_id(1),cell_id(2))

%% Pressure Animation
plot_animation(deltaCS,H,nz,cell_plot,time_plot_days,'Pressure_Data','\mathrm{Pressure~Head}','m',head_pressure,colors(1,:),cell_id(1),cell_id(2))

%% Moisture Animation
plot_animation(deltaCS,H,nz,cell_plot,time_plot_days,'Moisture_Data','\theta','cm^3\cdot cm^{-3}',waterContent,colors(2,:),cell_id(1),cell_id(2))

%% Pressure Profiles
profile_plots(6,time_plot_days,cell_plot(1),cell_plot(2),cell_id(1),cell_id(2),head_pressure,(cumsum(deltaCS)),'pressure_profiles','h','m');

%% Moisture Profiles
profile_plots(6,time_plot_days,cell_plot(1),cell_plot(2),cell_id(1),cell_id(2),waterContent,(cumsum(deltaCS)),'moisture_profiles','\theta','cm^3 \cdot cm^{-3}');

%% GW
[col_out,row_out] = find(lsurf == min(min(lsurf)),1,'first');
line_plot(time_plot_days,squeeze(GW(1,col_out,row_out,:)),'{h}_{\mathrm{GW}}','m',-TopBC*3600*1000,ChangeT/86400,'i','\mathrm{mm\cdot h^{-1}}','GW Table at Chosen Cell')
exportgraphics(gcf,'Output/GW_Table_at_Outlet.png','Resolution','600')

%% Hydrograph
line_plot(time_plot_days,-discharge*1000*3600,'Q','\mathrm{mm\cdot h^{-1}}',-TopBC*3600*1000,ChangeT/86400,'i','\mathrm{mm\cdot h^{-1}}','Hydrograph')
exportgraphics(gcf,'Output/GW_Hydrograph.png','Resolution','600')
% line_plot(time_plot_days,1000*discharge*(nx*dx*ny*dy),'Q','\mathrm{L\cdot s^{-1}}')

%% Storage-Discharge
scatter_plot(1000*GW_Storage/(dx*dy*numel(X)),-discharge*1000*3600,time_plot_days,'\bar{h}_{\mathrm{GW}}','Q','mm','mm \cdot h^{-1}','GW Storage-Discharge')
exportgraphics(gcf,'Output/GW_Storage_Discharge.png','Resolution','600')

%% Pressure Head
depths_field = [0.1 0.5 0.75];
field_plots(H,time_plot_days,'Head Pressure','P/{\gamma_w}','m',head_pressure,zone_coordinates(:,1),zone_coordinates(:,2),X,Y,depths_field,deltaCS,1,[0,90])

%% Level Plots
nx_values = [1  12  24];
ny_values = [1  30  60];
level_plots(time_plot_days,'Water Level Plots',GW,z0,nx_values,ny_values)

%% Potential Field
depths_field = [0.1 0.5 0.75];
field_plots(max(max(max(max(potential_head)))),time_plot_days,'Potential Head','P/{\gamma_w} + z','m',potential_head,zone_coordinates(:,1),zone_coordinates(:,2),X,Y,depths_field,deltaCS,0,[0,90])

%% Field Groundwater
field_plots(H,time_plot_days,'Groundwater_Field','GW','m',GW,Y(:),X(:),X,Y,1,1,0,[0,90])

%% 3D Field Groundwater Depth
field_plots(H,time_plot_days,'Groundwater_Field_3D','GW','m',GW,Y(:),X(:),X,Y,1,1,0,[45,60])

%% GW
intervals = 5;
rows = 2;
for i = 1:(intervals+1)
    subplot(rows,ceil((intervals+1)/rows),i);
    dx_time = length(time_plot_days)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    surf_plot(H,time_plot_days(idx_time),'h_{GW}','m',GW_original(:,:,idx_time),1,1,256,0.9,0,[0,90]);
end
% %% 3D Field Groundwater Table
% field_plots(max(max(max(max(GW_table)))),time_plot_days,'Groundwater_Table_3D','GW','m',GW_table,Y(:),X(:),X,Y,1,1,0,[45,60])
exportgraphics(gcf,'Output/GW_Depths_Surfs.png','Resolution','600')
%% Animation Plots
function plot_animation(deltaCS,H,nz,cell_plot,time_plot,name_plot,symbol,units,variable,color,y_cell,x_cell)
soil_depth = (cumsum(deltaCS));
% Video
figure(1)
obj = VideoWriter(['Output/Videos' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
for n=1:1:(size(variable,4))
    clf
    if n == 1
        t = 0;
    else
        t = time_plot(n); % days
    end
    vector_data = variable(:,cell_plot(1),cell_plot(2),n);
    x = linspace(max(max(max(max(vector_data)))),min(min(min(min(vector_data)))),nz);
    y = 0*soil_depth;
    plot(x,y,'LineWidth',4,'LineStyle','-','Color','k')
    %     plot(0*x,soil_depth,'LineWidth',4,'LineStyle','-','Color','k')
    hold on
    plot(vector_data,soil_depth,'k','LineWidth',2,'LineStyle','-','Color',color);
    % Marker Properties
    % 'marker','v','MarkerEdgeColor','black','MarkerSize',1,'MarkerFaceColor','red');
    xlabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex');
    ylabel('$z$ [m] - From the bottom','Interpreter','latex');
    axis tight
    set(gca,'TickDir','out');
    set(gca,'FontName','Garammond')
    set(gca,'FontSize',14);
    set(gca,'FontWeight','Bold','LineWidth', 1.5);
    set(gca,'TickLength',[0.02 0.01])
    xlim([min(min(min(min(variable)))),max(max(max(max(variable))))])
    grid off
    box on
    title(['t = ',num2str(t),' [days]   -   ycell [',num2str(y_cell),'], ','xcell [',num2str(x_cell),']'])
    f = getframe(gcf);
    writeVideo(obj,f);
end
drawnow
obj.close();
end


%% Field Plots
function field_plots(H,time_plot,name_plot,symbol,units,variable,yi,xi,X,Y,depths_field,deltaCS,plot_title,angle)
% Video
close all
figure(2)
clf
obj = VideoWriter(['Output/Videos' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 5;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
for n=1:1:(size(variable,4))
    if n == 1
        t = 0;
    else
        t = time_plot(n); % days
    end
    sgtitle(['t = ',num2str(t),' [days]'],'interpreter','latex','FontSize',16,'FontWeight','bold')

    for k = 1:length(depths_field)

        subplot(1,length(depths_field),k)
        dz_idx = find(depths_field(k) <= cumsum(deltaCS),1,'first');
        var = squeeze(variable(dz_idx,:,:,n));
        % If number of data equals the meshgrid size
        if numel(var) == numel(X)
            interp = var;
        else
            interp = idw_function([xi,yi],var(:),X,Y,1,inf,1);

        end
        vector_data = interp;
        surf_plot(H,time_plot,symbol,units,vector_data,1,0,64,0.9,0,angle);

        if plot_title == 1
            title(['$\Delta z$ = ',num2str(depths_field(k)),' [m]'],'Interpreter','latex')
        end
        f = getframe(gcf);
        writeVideo(obj,f);
    end
end
drawnow
obj.close();
end



%% Level Plots
function level_plots(time_plot,name_plot,GW,ground_el,nx_values,ny_values)
% Video
figure(4)
clf
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
obj = VideoWriter(['Output/Videos' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
for n=1:1:(size(GW,4))
    matrix = squeeze(GW(1,:,:,n));
    level_plot(time_plot(n),matrix,ground_el,nx_values,ny_values)
    f = getframe(gcf);
    writeVideo(obj,f);
end
end

%% Discharge Plot
function line_plot(time_plot,y_values,symbol,units,topBC,TimeBC,symbol_T,units_T,title_name)
close all
c = linspecer(10);

set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
plot(time_plot,y_values,'Color',c(2,:),'linewidth',2,'Marker','*')
xlabel('Elapsed Time [days]','Interpreter','latex')
ylabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
hold on


if (sum(topBC) ~= 0) && (sum(topBC(topBC > 0) + sum(topBC(topBC < 0)))) ~= 0
yyaxis right
set(gca,'ydir','reverse');
set(gca,'YColor','black')
bar(TimeBC(1:length(topBC)),topBC,'FaceColor',[0 128 128]/256,'EdgeColor',[.3 .3 .3],'LineWidth',1.2)
ylabel(['$ ' symbol_T '~[ \mathrm{ ' units_T '}]$'],'Interpreter','latex');
try
ylim([0,10*max(topBC)])
catch 
ylim([0,1])    
end
Flow = sum(y_values*(24*time_plot(1)));
Rain = sum(topBC*(24*TimeBC(2)));
end
title(title_name,'Interpreter','latex','FontSize',14);

end

%% Scatter Plot
function scatter_plot(variable_1,variable_2,time_plot_days,symbol_1,symbol_2,units_1,units_2,title_name)
close all
sz = 5;
c = linspace(1,length(variable_1),length(variable_1));
set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
scatter(variable_1,variable_2,sz,c,'filled');
colormap(['winter(',num2str(length(variable_1)),')'])
xlabel(['$ ' symbol_1 '~[ \mathrm{ ' units_1 '}]$'],'Interpreter','latex');
ylabel(['$ ' symbol_2 '~[ \mathrm{ ' units_2 '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
time_factor = 1/length(variable_1)*time_plot_days(end);
ticks = [1 length(variable_1)/5 2*length(variable_1)/5 3*length(variable_1)/5 4*length(variable_1)/5 5*length(variable_1)/5];
h = colorbar;
h.Ticks = (ticks);
h.TickLabels = [round(time_factor*ticks,2)];
ylabel(h,'Elapsed Time [days]','Interpreter','latex','FontSize',14);
h.TickDirection = 'out';
h.FontSize = 14;
hold on
box on
title(title_name,'Interpreter','latex','FontSize',14);
end
%% Elevation Video
function elevation_surf(time_plot,name_plot,symbol,units,variable)
% Video
close all
figure(2)
clf
obj = VideoWriter(['Output/Videos' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
nmax = 20;
for n=1:1:nmax
    sgtitle(['Bedrock Elevation [m]'],'interpreter','latex','FontSize',16,'FontWeight','bold')
    vector_data = variable;
    surf_plot(max(max(vector_data)),time_plot(1),symbol,units,vector_data,1,0,64,0.9,0,[0,90]);
    view(-(n-1)*360/nmax,(n-1)*90/nmax);
    f = getframe(gcf);
    writeVideo(obj,f);
end
drawnow
obj.close();
end

%% Profile Plots
function [handle] = profile_plots(intervals,time_plot,y_cell,x_cell,y_id,x_id,four_dim_data,z_cells,plot_name,symbol,units)
close all
for i = 1:(intervals+1)
    dx_time = length(time_plot)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    y_values = squeeze(four_dim_data(:,y_cell,x_cell,idx_time));
    c = linspecer(intervals+1);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
    handle = plot(y_values,z_cells,'Color',c(i,:),'linewidth',2);
    xlabel(['$',symbol,'~[\mathrm{',units,'}]$'],'Interpreter','latex')
    ylabel('Elevation from bedrock [m]','Interpreter','latex');
    set(gca,'Fontsize',14)
    set(gca,'FontName','garamond');
    set(gca,'TickDir','out');
    day_stg = 'd';
    label{i} = ['$ t = $',num2str(round(time_plot(idx_time),2)),' [',day_stg,']'];
    hold on
end
title(['y(cell) = [',num2str(y_id),'], ','x(cell) = [',num2str(x_id),']'])

legend(label,'Interpreter','latex','FontSize',12);
% %% 3D Field Groundwater Table
% field_plots(max(max(max(max(GW_table)))),time_plot_days,'Groundwater_Table_3D','GW','m',GW_table,Y(:),X(:),X,Y,1,1,0,[45,60])
exportgraphics(gcf,['Output/',plot_name,'.png'],'Resolution','600')
end