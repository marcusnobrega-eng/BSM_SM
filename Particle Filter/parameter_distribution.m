function parameter_distribution(parameters,correct_parameters,time_min,symbol_1,symbol_2,units_1,units_2)
% Developer: Marcus Nobrega, Ph.D
% variable_1: vector 
% time_plot: time vector in days associaed with var 1 and var 2
% symbol_1: symbol of var 1 (e.g., 'Q')
% symbol_2: symbol of var 2
% units_1: units of var 1 (e.g., 'mm/h')
% units_2: units of var 2
% title_name: name to plot at the top of the chart

color1 = [76, 73, 68]/255;
color2 = [255, 29, 147]/255;
sz = 12;
x_parameters = 1:1:length(correct_parameters);
plot(x_parameters,correct_parameters,'color',color1,'LineWidth',3,'LineStyle','--');
hold on
scatter(x_parameters,parameters,sz,color2);
xlabel(['$ ' symbol_1 '~[ \mathrm{ ' units_1 '}]$'],'Interpreter','latex');
ylabel(['$ ' symbol_2 '~[ \mathrm{ ' units_2 '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
legend('Correct','Particle Filter','interpreter','latex');
hold on
% scatter(0,observed_state,150,color2,'filled')
box on
title(strcat('t =~',num2str(round(time_min,2)),'~min'),'Interpreter','latex','FontSize',14);
% pause(0.001)
end