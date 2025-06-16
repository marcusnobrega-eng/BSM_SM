function level_plots(time_plot,name_plot,variable_1,var1_name,variable_2,var2_name,nx_values,ny_values,X,Y,Hmin,Hmax,y_label,y_label_units)
% Creates GW Table profiles
% Developer: Marcus Nobrega, Ph.D
% time_plot: time vector in days
% variable: 3D variable with third dimension being time
% name_plot: name of the file to be saved
% ground_el: bedrock elevation matrix
% nx_values: cells in x (cols) e.g., [1 4 nx], nx > 4
% ny_values: cells in y (rows) e.g., [2 3 ny], ny > 3
% X: x meshgrid [m]
% Y: y meshgrid [m]
% Hmin: min value of the plot
% Hmax: max value of the plot

% Video
figure(4)
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end
clf
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
obj = VideoWriter(['Output/Animations/' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
for n=1:1:(size(variable_1,3))
    matrix = squeeze(variable_1(:,:,n));
    level_plot(time_plot(n),matrix,variable_2,nx_values,ny_values,X,Y,var1_name,var2_name,Hmin,Hmax,y_label,y_label_units)
    f = getframe(gcf);
    writeVideo(obj,f);
end
end