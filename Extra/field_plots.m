% Field Plots
function field_plots(H,time_plot,name_plot,symbol,units,variable,angle,variable_name,X,Y)
% Developer: Marcus Nobrega, Ph.D
% Goal: Surfplot Videos of a 3D Array
% H: maximum value of the variable
% time_plot: vector with time in days
% name_plot: name of the plot to save
% symbol: e.g., h. You do not need to enter $ for instance
% units: e.g., m. You do not need to enter [m].
% variable: a 3-D variable
% angle: [0,90] for a plan view
% X: meshgrid for X values in [m]
% Y: meshgrid for Y values in [m]

close all
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end
figure(2)
clf
obj = VideoWriter(['Output/Animations' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 5;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
for n=1:1:(size(variable,3))
    if n == 1
        t = 0;
    else
        t = time_plot(n); % days
    end
    sgtitle(['t = ',num2str(t),' [days]'],'interpreter','latex','FontSize',16,'FontWeight','bold')

    for k = 1:1

        subplot(1,1,k)
        var = squeeze(variable(:,:,n));
        % If number of data equals the meshgrid size
        vector_data = var;
        surf_plot(H,time_plot,symbol,units,vector_data,1,0,64,0.9,0,angle,X,Y);

        title([variable_name],'Interpreter','latex')
        f = getframe(gcf);
        writeVideo(obj,f);
    end
end
drawnow
obj.close();
end
