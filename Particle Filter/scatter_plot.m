%% Scatter Plot
function scatter_plot(variable_1,variable_2,time_plot_days,symbol_1,symbol_2,units_1,units_2,title_name)
% Developer: Marcus Nobrega, Ph.D
% variable_1: vector 
% varialbe_2: vector
% time_plot: time vector in days associaed with var 1 and var 2
% symbol_1: symbol of var 1 (e.g., 'Q')
% symbol_2: symbol of var 2
% units_1: units of var 1 (e.g., 'mm/h')
% units_2: units of var 2
% title_name: name to plot at the top of the chart

close all
sz = 5;
c = linspace(1,length(variable_2),length(variable_2));
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