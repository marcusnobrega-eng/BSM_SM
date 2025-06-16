%% Scatter Plot
function particle_distribution(estimated_state,probability,observed_state,time_min,symbol_1,symbol_2,units_1,units_2,legend_1,legend_2,ylimval)
% Developer: Marcus Nobrega, Ph.D
% variable_1: vector 
% time_plot: time vector in days associaed with var 1 and var 2
% symbol_1: symbol of var 1 (e.g., 'Q')
% symbol_2: symbol of var 2
% units_1: units of var 1 (e.g., 'mm/h')
% units_2: units of var 2
% legend_1: legend of first entry
% legend_2: legend of second entry

color1 = [76, 73, 68]/255;
color2 = [255, 29, 147]/255;
sz = 24;
% if time_min == 0
% set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, .25, .5, .5]);
% else
%     clf
% end
c = linspace(1,length(probability),length(probability));

scatter(probability,estimated_state,sz,c);
colormap(['turbo(',num2str(length(estimated_state)),')'])
xlabel(['$ ' symbol_1 '~[ \mathrm{ ' units_1 '}]$'],'Interpreter','latex');
ylabel(['$ ' symbol_2 '~[ \mathrm{ ' units_2 '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
hold on
scatter(0,observed_state,150,color2,'filled')
box on
h = colorbar;
ylabel(h,'Particle Index','Interpreter','latex','FontSize',12);
h.TickDirection = 'out';
h.FontSize = 12;
title(strcat('t =~',num2str(round(time_min,2)),'~min'),'Interpreter','latex','FontSize',14);
legend(legend_1,legend_2,'interpreter','latex');
ylim(ylimval);
% if max(estimated_state) == min(estimated_state)
%     ylim([(min(estimated_state)), (max(estimated_state)+0.5*(min(estimated_state))+0.1)])
% else
%     min_val = min(min(estimated_state),observed_state);
%     max_val = max(min(estimated_state),observed_state);
%     ylim([min_val, max_val]);
% end

% % Adding labels to the particles
% b = string(1:1:length(probability));
% x = probability;
% y = estimated_state;
% c = cellstr(b);
% dx = 0.005; dy = 0.000; % displacement so the text does not overlay the data points
% text(x+dx, y+dy, c,"FontSize",8,'Interpreter','latex');

% pause(0.001)
end