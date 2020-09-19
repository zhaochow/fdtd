function plotE2(evo_time)
t = 0;

vName = input('Enter video file name>', 's');
video_export = VideoWriter(char(vName+".avi"));
video_export.FrameRate = 60;
video_export.Quality = 100;
open(video_export)
input('Press any key to start>', 's')

[grid_sizeX, grid_sizeY, t_final] = size(evo_time);
Xaxis = 1:grid_sizeX;

grid_sizeY = 2 * 11;

figure
h = plot(Xaxis,evo_time(:,round(grid_sizeY/2),t+1));
% h = plot(Xaxis,evo_time(round(grid_sizeY/2),:,t+1));
% h = surf(evo_time(:,:,t+1));
ylim([-1 1]);
xlabel('Y step');
ylabel('E [V/m]');
title('Electric field at X = 370')

while t < t_final-1
    t = t+1;
    
    set(h, 'YData', evo_time(:,round(grid_sizeY/2),t+1));
%     set(h, 'YData', evo_time(round(grid_sizeY/2),:,t+1));
%     set(h, 'ZData', evo_time(:,:,t+1));
    drawnow
    M(t) = getframe(gcf);  %Video frames
    pause(0.01)
end

open(video_export)
writeVideo(video_export, M);
close(video_export);

end