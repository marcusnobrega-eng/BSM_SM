function line_plot(x_value,x_symbol,x_units,y_values,symbol,units,topBC,TimeBC,symbol_T,units_T,title_name)
% Create a hydrograph in left axis with a boundary condition in the reverse right axis
% x_value: time vector 
% x_symbol: symbol of x e.g., 'x'
% y_values: discharge or other vector
% symbol: symbol of the y_values (e.g., 'Q')
% units: units of the y_values (e.g., '\mathrm{mm \cdot h^{-1}}
% topBC: rainfall or other top B.C vector to plot with the y_values
% symbol_T: symbol of the top B.C. (e.g., 'i')
% units_T: units of the top B.C.
% title_name: Name to plot at the top of the graph

close all
c = linspecer(10);

matrix = topBC;
topBC = 0;
% Adjusting TopBC
if size(matrix,3) > 1
    for i = 1:size(matrix,3)
        data = matrix(:,:,i);
        topBC(1,i) = mean(data(:));
    end
else % Single value per area
    topBC = matrix;
end

set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
plot(x_value,y_values,'Color',c(1,:),'linewidth',2)
xlabel(['$ ' x_symbol '~[ \mathrm{ ' x_units '}]$'],'Interpreter','latex')
ylabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
hold on
try
    ylim([0,1.5*max(y_values)])
catch
    ylim([min(y_values),1.5*max(y_values)])
end


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
end

title(title_name,'Interpreter','latex','FontSize',14);

end