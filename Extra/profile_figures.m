function profile_figures(intervals,time_plot,matrix_3D,X,Y,nx_values,ny_values,plot_name,symbol,units)
% Developer: Marcus Nobrega, Ph.D.
% Goal: Plot Profile plots of a 3D variable
% intervals: number of finite intervals of time
% time_plot: time in days associated with the data
% matrix_3D: 3D matrix with the data
% nx_values: indexes of cells which data will be ploted
% ny_values: indexes of cells which data will be ploted
% plot_name: name of the plot
% symbol: symbol of the variable
% units: units of the variable

m = 0;
n_plots = max(length(nx_values),length(ny_values))*2;
close all
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
for j = 1:length(nx_values)
    m = m + 1;
    subplot(ceil(n_plots/2),2,(m-1)*2+1);
    for i = 1:(intervals+1)
        dx_time = length(time_plot)/intervals;
        idx_time = floor((i-1)*dx_time);
        if i == 1
            idx_time = 1;
        end
        y_values = squeeze(matrix_3D(:,nx_values(j),idx_time));
        c = linspecer(intervals+1);
        x_values = Y(:,nx_values(j));
        set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
        if mod(i,2) == 0
            handle = plot(x_values,y_values,'Color',c(i,:),'linewidth',2,'LineStyle','-');
        else
            handle = plot(x_values,y_values,'Color',c(i,:),'linewidth',2,'LineStyle',':');
        end
        ylabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex')
        xlabel('Distance alongside $y$ [m]','interpreter','latex')
        set(gca,'Fontsize',14)
        set(gca,'FontName','garamond');
        set(gca,'TickDir','out');
        day_stg = 'd';
        hold on
        label{i} = ['$ t = $',num2str(round(time_plot(idx_time),2)),' [',day_stg,']'];
        if m == 1
            legend(label,'interpreter','latex','Location','bestoutside');
        end
    end
    title_lbel = ['$ x = ',num2str(round(nx_values(j),2)),'~\mathrm{m}$'];
    title(title_lbel,'interpreter','latex');
end

m = 0;
for j = 1:length(ny_values)
    m = m + 1;
    subplot(ceil(n_plots/2),2,(m)*2);
    for i = 1:(intervals+1)
        dx_time = length(time_plot)/intervals;
        idx_time = floor((i-1)*dx_time);
        if i == 1
            idx_time = 1;
        end
        y_values = squeeze(matrix_3D(ny_values(j),:,idx_time));
        c = linspecer(intervals+1);
        x_values = X(ny_values(j),:);
        set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
        if mod(i,2) == 0
            handle = plot(x_values,y_values,'Color',c(i,:),'linewidth',2,'LineStyle','-');
        else
            handle = plot(x_values,y_values,'Color',c(i,:),'linewidth',2,'LineStyle',':');
        end
        ylabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex')
        xlabel('Distance alongside $x$ [m]','interpreter','latex')
        set(gca,'Fontsize',14)
        set(gca,'FontName','garamond');
        set(gca,'TickDir','out');
        day_stg = 'd';
        hold on
    end
    title_lbel = ['$ y = ',num2str(round(ny_values(j),2)),'~\mathrm{m}$'];
    title(title_lbel,'interpreter','latex');    
end
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
exportgraphics(gcf,['Output/Figures/',plot_name,'.png'],'Resolution',600)
end

