function [handle] = level_plot(time_plot,variable_1,variable_2,nx_values,ny_values,X,Y,var1_name,var2_name,Hmin,Hmax,y_label,y_label_units)
% Creates GW Table profiles
% Developer: Marcus Nobrega, Ph.D
% time_plot: time vector in days
% variable_1: 2D matrix with values per cell
% name_plot: name of the file to be saved
% variable_2: typycally the bedrock elevation matrix
% nx_values: cells in x (cols) e.g., [1 4 nx], nx > 4
% ny_values: cells in y (rows) e.g., [2 3 ny], ny > 3
% X: x meshgrid [m]
% Y: y meshgrid [m]
% var1_name: name of variable 1 
% var2_name: name of variable 2 
% Hmin: minimum value of the plots
% Hmax: maximum value of the plots
% y_label: e.g., 'head'
% y_label_units: e.g. [m]
time_frame = time_plot(1); % days
m = 0;
clf
n_plots = length(ny_values) + length(nx_values);
sgtitle(['t = ',num2str(time_frame),' [days]'],'interpreter','latex','FontSize',18,'FontWeight','bold')
for i = 1:length(nx_values)
    m = m + 1;
    handle{m} = subplot(ceil(n_plots/2),2,(m-1)*2+1);
    g_el = variable_2(:,nx_values(i));
    h_gw = variable_1(:,nx_values(i));
    x_val = Y(:,nx_values(i));
    title_val = X(1,nx_values(i));
    plot(x_val,g_el,'Color','black','LineWidth',2);
    hold on
    plot(x_val,g_el + h_gw,'Color','blue','LineWidth',3,'LineStyle',':');
    xlabel('y [m]','FontSize',14,'Interpreter','latex')
    ylabel(['$ ' y_label '~[ \mathrm{ ' y_label_units '}]$'],'Interpreter','latex')
    title(['x = ',num2str(title_val),' [m]'],'interpreter','latex')
    set(gca,'TickDir','out');
    set(gca,'FontName','Garammond')
    set(gca,'FontSize',14);
    grid off
    box on
    axis tight
    try
%         ylim([min(min(variable_2)),max(max(variable_1 + variable_2))])
        ylim([Hmin,Hmax])     
    catch ME

    end
    hold on
    if m == 1
        legend(['$',var2_name,'$'],['$',var1_name,'$'],'interpreter','latex')
    end
end
m = 0;
for i = 1:length(ny_values)
    m = m + 1;
    handle{m} = subplot(ceil(n_plots/2),2,(m)*2);
    g_el = variable_2(ny_values(i),:);
    h_gw = variable_1(ny_values(i),:);
    y_val = X(ny_values(i),:);
    title_val = Y(ny_values(i),1);
    plot(y_val,g_el,'Color','black','LineWidth',2);
    hold on
    plot(y_val,g_el + h_gw,'Color','blue','LineWidth',3,'LineStyle',':');
    ylabel(['$ ' y_label '~[ \mathrm{ ' y_label_units '}]$'],'Interpreter','latex')    
    title(['y = ',num2str(title_val),' [m]'],'interpreter','latex')
    xlabel('x [m]','FontSize',14,'Interpreter','latex')    
    set(gca,'TickDir','out');
    set(gca,'FontName','Garammond')
    set(gca,'FontSize',14);
    grid off
    box on
    axis tight
    try
%         ylim([min(min(variable_2)),max(max(variable_1 + variable_2))])
        ylim([Hmin,Hmax])    
    catch ME

    end
    hold on  
end
hold on
drawnow
end
