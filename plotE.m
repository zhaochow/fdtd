function plotE(evo_time)
t = 0;
t_final = size(evo_time,3);

vName = input('Enter video file name>', 's');
video_export = VideoWriter(char(vName+".avi"));
video_export.FrameRate = 60;
video_export.Quality = 100;
open(video_export)
input('Press any key to start>', 's')

% figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
figure
h = imagesc(db(abs(evo_time(:,:,t+1))), [-50 0]);
metal_coord = [500/2 - 50-0.5,500/2-10-2.5,...
    2*50+1,3];
rectangle('Position',metal_coord,'LineWidth',1,'LineStyle','-', ...
    'EdgeColor', 'w','FaceColor', 'w')
hold on
colormap(jet)
axis image

while t < t_final
    t = t+1;
    
    set(h, 'CData', db(abs(evo_time(:,:,t))));
    drawnow
    M(t) = getframe;  %Video frames
    pause(0.01)
end

open(video_export)
writeVideo(video_export, M);
close(video_export);


end