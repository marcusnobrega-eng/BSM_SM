%% Animation Plots
function plot_animation(deltaCS,nz,cell_plot,time_plot,name_plot,symbol,units,variable,color,y_cell,x_cell)
% Developer: Marcus Nobrega, Ph.D.
% Goal: Generate videos of plots of states for UZ 
% deltaCS: vector with the distance of each UZ node from the bedrock
% nz: number of UZ nodes
% cell_plot: index of the cell to plot [<= ny, <= nx]
% time_plot: time vector in days
% name_plot: name of the output file
% symbol: symbol of the variable plotted
% units: units of the variable plotted
% variable: 4-D array with the state
% color: color of the line 
% y_cell: y coordinate of the cell
% x_cell: x coordinate of the cell

soil_depth = (cumsum(deltaCS));
% Video
figure(1)
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end

obj = VideoWriter(['Output/Animations/' name_plot '.avi'],'Motion JPEG AVI');
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
    title(['t = ',num2str(round(t,2)),' [days]   -   y = ',num2str(round(y_cell,2)),', m ','x = ',num2str(round(x_cell,2)),' m'])
    f = getframe(gcf);
    writeVideo(obj,f);
end
drawnow
obj.close();
end