clear; 
close all;

%% Some useful variables

grid_sizeX = 500; % Width of the simulation grid
grid_sizeY = 500; % Height of the simulation grid

fields = zeros(grid_sizeY,grid_sizeX,3); % Ez(l,m,n),Hx(l,m+1/2,n),Hy(l+1/2,m,n)
EZ = 1;
HX = 2;
HY = 3;
param = zeros(grid_sizeY,grid_sizeX,3); % deltaT/(eps*deltaX),deltaT/(mu*deltaX)
EPS = 1;
MU = 2;


%% Parameters

c = 299792458; % Speed of light
eps_vac = 8.85418782e-12; % Permittivity of the vacuum
mu_vac = 4*pi*1e-7; % Permeability of the vacuum
ref_idx_diamond = 2.417; % Refraction index of the diamond
eps_diamond = (ref_idx_diamond^2)*eps_vac;     %Only valid in THz region (10^12)
% freq = 500e12; % Frequency of the excitation
freq = 2.45e9;
lambda = c/freq; % Corresponding wavelength
sourceX = grid_sizeX/2;
sourceY = grid_sizeY/2;

t = 0;
t_final = 1200;
% evo_time = zeros(grid_sizeY,grid_sizeX,t_final+1); % Variable for saving the simulation results

% deltaX = lambda/50;    % Boosted precision
deltaX = lambda/10;
deltaT = deltaX/(c*sqrt(2));

c_mu = deltaT/(mu_vac*deltaX);
c_eps = deltaT/(eps_vac*deltaX);
c_eps_conductor = deltaT/(1000000000000*deltaX);

c_eps_diamond = deltaT/(eps_diamond*deltaX);

% Initialize parameters
param(:,:,EPS) = c_eps;
param(:,:,MU) = c_mu;

grid_unitY = grid_sizeY/10;
grid_unitX = grid_sizeX/10;


%% Choose one of the following simulations
% WAVEGUIDE

% guide_width = ceil(0.2/deltaX);
% 
% 
% y1 = round(grid_sizeY/2 - (guide_width)/2); % Upper bound
% y2 = round(grid_sizeY/2 + (guide_width)/2); % Lower bound
% 
% for j = 1:grid_sizeX
%     param(y1,j,EPS) = c_eps_conductor;
%     param(y2,j,EPS) = c_eps_conductor;
% end

% DIFFRACTION

% x1 = round(5*grid_sizeX/8);
% 
% for i = 4*grid_unitY:grid_sizeY-4*grid_unitY
%     param(i,x1,EPS) = c_eps_conductor;
% end

% DOUBLE-SLIT

% x1 = round(5*grid_sizeX/9);
% 
% for i = 1:grid_sizeY
%     param(i,x1,EPS) = c_eps_conductor;
% end
% param(grid_sizeY/2-30:grid_sizeY/2-28,x1,EPS) = c_eps;
% param(grid_sizeY/2+28:grid_sizeY/2+30,x1,EPS) = c_eps;

% PERFECT CONDUCTOR BLOCK

% for i = grid_sizeX/2 - grid_unitX:grid_sizeX/2 + grid_unitX
%     for j = grid_sizeY/2 - grid_unitX:grid_sizeY/2 + grid_unitX
%          param(i,j,EPS) = c_eps_conductor;
%     end
% end

% TRANSMISSION - DIAMOND

% rect_coord = [grid_sizeX/2 - grid_unitX/2,grid_sizeY/2 - 2*grid_unitY,grid_unitX/2,4*grid_unitY]; % Coordinates of the diamond rectangle
% for i = grid_sizeX/2 - grid_unitX/2:grid_sizeX/2
%     for j = grid_sizeY/2 - 2*grid_unitY:grid_sizeY/2 + 2*grid_unitY
%          param(j,i,EPS) = c_eps_diamond;
%     end
% end

% TROPOSPHERIC REFRACTION

% step = 100;
% h = (grid_sizeY-1)*step:-step:0;
% H = 7.35e3;
% Nt = 315*exp(-h/H);
% for i = 1:grid_sizeY
%     ref_idx_tropo = 1 + 1e-6*Nt(i); % Refraction index
%     eps_tropo = (ref_idx_tropo^2)*eps_vac;
%     for j = 1:grid_sizeX
%          param(i,j,EPS) = c_eps_diamond;
%     end
% end

% PARABOLA REFLECTOR

% focus = 40;
% maxX = round(sqrt(round(2*focus)*4*focus)); % parabola ends 1 focus further
% x = -maxX:maxX;
% y = x.^2/(4*focus);
% 
% reflectorX = sourceX + x;
% reflectorY = sourceY + focus - round(y);
% 
% for i = 1:length(reflectorX)-1
%     y1 = reflectorY(i);
%     x1 = reflectorX(i);
%     y2 = reflectorY(i+1);
%     x2 = reflectorX(i+1);
%     param(y1,x1,EPS) = c_eps_conductor;
%     param(y2,x2,EPS) = c_eps_conductor;
%     if y2 - y1 == -1
%         param(y1,x2,EPS) = c_eps_conductor;
%     elseif y2 - y1 == 1
%         param(y2,x1,EPS) = c_eps_conductor;
%     elseif y2 - y1 == -2
%         param(y1-1,x2,EPS) = c_eps_conductor;
%     elseif y2 - y1 == 2
%         param(y1+1,x1,EPS) = c_eps_conductor;
%     end
% end

%% Plotting and recording

figure()

% Plot 2D figure with dB scale (min,max) = (-50,0)

h = imagesc(db(abs(fields(:,:,EZ))), [-50 0]);
hold on
% rectangle('Position',rect_coord,'LineWidth',2,'LineStyle','--', 'EdgeColor', 'w')
colormap(jet)
axis image % Same x:y aspect ratio

focus = 40;
maxX = round(sqrt(round(2*focus)*4*focus)); % parabola ends 1 focus further
x = -maxX:maxX;
y = x.^2/(4*focus);

reflectorX = sourceX + x;
reflectorY = sourceY + focus - round(y);

for i = 1:length(reflectorX)-1
    y1 = reflectorY(i);
    x1 = reflectorX(i);
    y2 = reflectorY(i+1);
    x2 = reflectorX(i+1);
    param(y1,x1,EPS) = c_eps_conductor;
    rectangle('Position',[x1-0.5,y1-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    param(y2,x2,EPS) = c_eps_conductor;
    rectangle('Position',[x2-0.5,y2-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    if y2 - y1 == -1
        param(y1,x2,EPS) = c_eps_conductor;
        rectangle('Position',[x2-0.5,y1-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    elseif y2 - y1 == 1
        param(y2,x1,EPS) = c_eps_conductor;
        rectangle('Position',[x1-0.5,y2-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    elseif y2 - y1 == -2
        param(y1-1,x2,EPS) = c_eps_conductor;
        rectangle('Position',[x2-0.5,y1-1-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    elseif y2 - y1 == 2
        param(y1+1,x1,EPS) = c_eps_conductor;
        rectangle('Position',[x1-0.5,y1+1-0.5,1,1],'LineWidth',2,'LineStyle','-', 'EdgeColor', 'w', 'FaceColor', 'w')
    end
end

% Visualize in 3D

% h = surf(abs(fields(:,:,EZ)));
% set(h,'LineStyle','none')
% zlim([0 0.5])
% caxis([0 0.5])
% drawnow

% Record simulation

% vName = input('Enter video file name>', 's');
% video_export = VideoWriter(char(vName+".avi"));
% video_export.FrameRate = 60;
% video_export.Quality = 100;
% open(video_export)
% input('Press any key to start>', 's')

%% Loop

source_coeff = 2*pi*freq*deltaT; % Sinusoidal excitation

tic
while t < t_final
    matY = 2:grid_sizeY-1; % Do not take edges
    matX = 2:grid_sizeX-1; % Do not take edges    
    fields(matY,matX, HX) = fields(matY,matX, HX) ...
     - param(matY,matX,MU).*(fields(matY,matX+1, EZ) - fields(matY,matX,EZ)); % Update Hx
    fields(matY,matX, HY) = fields(matY,matX, HY) ...
     + param(matY,matX,MU).*(fields(matY+1,matX, EZ) - fields(matY,matX,EZ)); % Update Hy
    fields(matY,matX, EZ) = fields(matY,matX, EZ) ...
     + param(matY,matX,EPS).*(fields(matY,matX,HY) - fields(matY-1,matX,HY) - fields(matY,matX,HX) + fields(matY,matX-1,HX)); % Update Ez
    
%     fields(sourceY,sourceX,1) = exp(-(t-50)^2/100);
    if(t*deltaT <= 1/freq) % Time limited sinusoidal excitation
        fields(sourceY,sourceX,1) = sin(source_coeff*t);
    elseif (t*deltaT <= 3/freq)
        fields(sourceY,sourceX,1) = 0;
    end
%     fields(sourceY,sourceX,EZ) = sin(source_coeff*t); % Sinusoidal excitation
    
%     evo_time(:,:,t+1) = fields(:,:,EZ); % Save Ez values
    t = t+1; % Increment time
    
    % Draw 2D image
    
    set(h, 'CData', db(abs(fields(:,:,EZ))));
    drawnow
    
    % Draw 3D view
    
%     h.ZData = abs(fields(:,:,EZ));

    % Record simulation
    
%     M(t) = getframe;  %Video frames
end
toc

% Save video

% open(video_export)
% writeVideo(video_export, M);
% close(video_export);

%% Power visualization

% Point over time

% time = 0:t_final;
% m = grid_sizeX/2 - grid_unitX;
% n = grid_sizeY/2 - grid_unitX;
% figure
% plot(time,squeeze(evo_time(480,260,:)));

% Along X axis

% Xaxis = 1:grid_sizeX;
% figure
% plot(Xaxis,evo_time(grid_sizeY/2,:,71));

% Ez and power density along X axis
% Valid only for sinusoidal excitation, need multiple periods in the time window
% 
% Xaxis = 1:grid_sizeX;
% Ez = zeros(1,grid_sizeX);
% power = zeros(1,grid_sizeX);
% for i = 1:grid_sizeX
%     Ez(i) =
%     rms(squeeze(evo_time(grid_sizeY/2,i,t_final-20:t_final))*sqrt(2));
%     power(i) = Ez(i)^2/(2*120*pi); % Power density
% end
% 
% figure
% plot(Xaxis,Ez);
% title("Electric field")
% figure
% semilogy(Xaxis,power);
% title("Power density")

