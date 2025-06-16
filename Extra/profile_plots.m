%% Profile Plots
function [handle] = profile_plots(intervals,time_plot,y_cell,x_cell,y_id,x_id,four_dim_data,z_cells,plot_name,symbol,units)
% Developer: Marcus Nobrega, Ph.D.
% Goal: Plot Profile plots of a 4D variable
% intervals: number of finite intervals of time
% y_cell: cell which data is ploted (index)
% x_cell: cell which data is ploted (index)
% y_id: index of the 4D data
% x_id: index of the 4D data
% four_dim_data: 4D daa
% z_cells: uz height from the bottom
% plot_name: name of the plot
% symbol: symbol of the variable
% units: units of the variable

close all
try 
    mkdir Output/Figures
catch 
    warning('Folder (Figures) already created')
end
for i = 1:(intervals+1)
    dx_time = length(time_plot)/intervals;
    idx_time = floor((i-1)*dx_time);
    if i == 1
        idx_time = 1;
    end
    y_values = squeeze(four_dim_data(:,y_cell,x_cell,idx_time));
    c = linspecer(intervals+1);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
    if mod(i,2) == 0
        handle = plot(y_values,z_cells,'Color',c(i,:),'linewidth',2,'LineStyle','-');
    else
        handle = plot(y_values,z_cells,'Color',c(i,:),'linewidth',2,'LineStyle',':');
    end
    xlabel(['$',symbol,'~[\mathrm{',units,'}]$'],'Interpreter','latex')
    ylabel('Elevation from bedrock [m]','Interpreter','latex');
    set(gca,'Fontsize',14)
    set(gca,'FontName','garamond');
    set(gca,'TickDir','out');
    day_stg = 'd';
    label{i} = ['$ t = $',num2str(round(time_plot(idx_time),2)),' [',day_stg,']'];
    hold on
end
title(['y = ',num2str(y_id),' m, ','x = ',num2str(x_id),' m'])

legend(label,'Interpreter','latex','FontSize',12);
% %% 3D Field Groundwater Table
% field_plots(max(max(max(max(GW_table)))),time_plot_days,'Groundwater_Table_3D','GW','m',GW_table,Y(:),X(:),X,Y,1,1,0,[45,60])
exportgraphics(gcf,['Output/Figures/',plot_name,'.png'],'Resolution','600')
end