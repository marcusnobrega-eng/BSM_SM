%% Scatter Plot
function scatter_plot_scaled(variable_1,variable_2,symbol_1,symbol_2,units_1,units_2,scale_relationship,title_name)
% Developer: Marcus Nobrega, Ph.D
% variable_1: vector 
% varialbe_2: vector
% symbol_1: symbol of var 1 (e.g., 'Q')
% symbol_2: symbol of var 2
% units_1: units of var 1 (e.g., 'mm/h')
% units_2: units of var 2
% scale_relationship: vector with z values of variable_1 and variable_2
% title_name: name to plot at the top of the chart

% close all
szmax = 40;
szmin = 5;
sz = szmin + (scale_relationship - min(scale_relationship))*(szmax - szmin)/(max(scale_relationship) - min(scale_relationship));
c = scale_relationship;
% set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
scatter(variable_1,variable_2,sz,c,'filled');
colormap(['cool(',num2str(length(variable_1)),')'])
xlabel(['$ ' symbol_1 '~[ \mathrm{ ' units_1 '}]$'],'Interpreter','latex');
ylabel(['$ ' symbol_2 '~[ \mathrm{ ' units_2 '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
scale_factor = scale_relationship;
max_scale = max(scale_factor);
ticks = [0 max_scale/5 2*max_scale/5 3*max_scale/5 4*max_scale/5 5*max_scale/5];
h = colorbar;
h.Ticks = (ticks);
h.TickLabels = round(100*[ticks],2);
ylabel(h,'Probability [%]','Interpreter','latex','FontSize',14);
h.TickDirection = 'out';
h.FontSize = 14;
hold on
box on
title(title_name,'Interpreter','latex','FontSize',14);
end