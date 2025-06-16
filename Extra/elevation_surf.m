%% Elevation Video
function elevation_surf(time_plot,name_plot,symbol,units,variable,X,Y)
% Video
close all
try 
    mkdir Output/Animations
catch 
    warning('Folder (Animations) already created')
end
figure(2)
clf
obj = VideoWriter(['Output/Videos' name_plot '.avi'],'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
nmax = 20;
for n=1:1:nmax
    sgtitle(['Bedrock Elevation [m]'],'interpreter','latex','FontSize',16,'FontWeight','bold')
    vector_data = variable;
    surf_plot(max(max(vector_data)),time_plot(1),symbol,units,vector_data,1,0,64,0.9,0,[0,90],X,Y);
    view(-(n-1)*360/nmax,(n-1)*90/nmax);
    f = getframe(gcf);
    writeVideo(obj,f);
end
drawnow
obj.close();
end