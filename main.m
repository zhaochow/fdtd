clear; 
close all;

%% Parameters

t = 0;
t_final = 300;

%% Some useful variables

EZ = 1;
HX = 2;
HY = 3;

EPS = 1;
MU = 2;
SIG = 3;

% sim = get_source_sim();
sim = get_beamforming_sim();

evo_time = zeros(sim.grid_sizeY,sim.grid_sizeX,t_final+1);


%% Plotting and recording

figure()

h = imagesc(db(abs(sim.fields(:,:,EZ))), [-50 0]);
hold on
if(isfield(sim,'rectangles'))
    for i = 1:size(sim.rectangles,1)
        rectangle('Position',sim.rectangles(i,:),'LineWidth',2,'LineStyle','--', 'EdgeColor', 'w')
        hold on
    end
end
colormap(jet)
axis image


% vName = input('Enter video file name>', 's');
% video_export = VideoWriter(char(vName+".avi"));
% video_export.FrameRate = 60;
% video_export.Quality = 100;
% open(video_export)
% input('Press any key to start>', 's')

%% Loop

matY = 2:sim.grid_sizeY-1;
matX = 2:sim.grid_sizeX-1;
tic
while t < t_final
    
    sim.fields(matY,matX, HX) = sim.fields(matY,matX, HX) ...
     - sim.param(matY,matX,MU).*(sim.fields(matY,matX+1, EZ) - sim.fields(matY,matX,EZ));
    sim.fields(matY,matX, HY) = sim.fields(matY,matX, HY) ...
     + sim.param(matY,matX,MU).*(sim.fields(matY+1,matX, EZ) - sim.fields(matY,matX,EZ));
    sim.fields(matY,matX, EZ) = sim.fields(matY,matX, EZ) ...
     + sim.param(matY,matX,EPS).*(sim.fields(matY,matX,HY) - sim.fields(matY-1,matX,HY) - sim.fields(matY,matX,HX) + sim.fields(matY,matX-1,HX));
 
    for i = 1:size(sim.source_idx,1)
        sim.fields(sim.source_idx(i,1),sim.source_idx(i,2),EZ) = sin(sim.source_coeff(i)*t + sim.source_phase(i));
    end

    evo_time(:,:,t+1) = sim.fields(:,:,EZ);
    t = t+1;
    
    set(h, 'CData', db(abs(sim.fields(:,:,EZ))));
    drawnow
%     h.ZData = abs(fields(:,:,EZ));
%     M(t) = getframe;  %Video frames
%     pause(0.01) % tic-toc with imagesc is 13 ms meaning 77 Hz is max refresh rate
end
toc

% open(video_export)
% writeVideo(video_export, M);
% close(video_export);
%% Visualization

time = 3000:4500;

% tmp = zeros(1,length(time));
a1 = 50; b1 = 375;
a = a1-16; b = b1-15;
c = 1; d = 1;
tmp = []; tmp2 = [];
% filt = 1/5*ones(1,5);
% for i = 1:size(evo_time,1)
%     for j = 1:size(evo_time,2)
%         tmp2(i,j,:) = conv(squeeze(evo_time(i,j,:)),filt);
%     end
% end
for i = 1:length(time)
    if mod(i,8) == 0
        a = a + c;
        if a > a1+15
            a = a1+15;
            b = b + d;
            c = -c;
        elseif a < a1-15
            a = a1-15;
            b = b + d;
            c = -c;
        end
        if b > b1+15
            b = b1+15;
            d = -d;
        elseif b < b1-15
            b = b1-15;
            d = -d;
        end
%         tmp = [tmp tmp2(a,b,i)];
%         tmp(i) = mean(mean(evo_time(:,:,i)));
%         tmp = [tmp mean(evo_time(a,b,i-round(length(time)/30)+1:i))];
        tmp = [tmp rms(squeeze(evo_time(a,b,i-7:i)))];
%         tmp(i) = evo_time(a,b,i);
    end
%     tmp(i) = evo_time(a,b,1000);
end

figure
% plot(time,squeeze(db(abs(evo_time(460,385,:)))));
% o = abs(evo_time(450,375,:));
plot(squeeze(db(abs(tmp(1:1:end)))));
% viz_point(time, evo_time, X, Y);
% viz_lineX(evo_time);
% viz_power();
