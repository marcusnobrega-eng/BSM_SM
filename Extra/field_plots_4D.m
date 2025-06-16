%% Field Plots
function field_plots_4D(H,time_plot,name_plot,symbol,units,variable,yi,xi,X,Y,depths_field,deltaCS,plot_title,angle)
% Video
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
for n=1:1:(size(variable,4))
    if n == 1
        t = 0;
    else
        t = time_plot(n); % days
    end
    sgtitle(['t = ',num2str(t),' [days]'],'interpreter','latex','FontSize',16,'FontWeight','bold')

    for k = 1:length(depths_field)

        subplot(1,length(depths_field),k)
        dz_idx = find(depths_field(k) <= cumsum(deltaCS),1,'first');
        var = squeeze(variable(dz_idx,:,:,n));
        % If number of data equals the meshgrid size
        if numel(var) == numel(X)
            interp = var;
        else
            interp = idw_function([xi,yi],var(:),X,Y,1,inf,1);
        end
        vector_data = interp;
        surf_plot(H,time_plot,symbol,units,vector_data,1,0,64,0.9,0,angle,X,Y);

        if plot_title == 1
            title(['$\Delta z$ = ',num2str(depths_field(k)),' [m]'],'Interpreter','latex')
        end
        f = getframe(gcf);
        writeVideo(obj,f);
    end
end
drawnow
obj.close();
end