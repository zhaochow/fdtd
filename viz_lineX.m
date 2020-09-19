function viz_lineX(evo_time)
% Along X axis
[grid_sizeX, grid_sizeY, t_end] = size(evo_time);
Xaxis = 1:grid_sizeX;
figure
plot(Xaxis,evo_time(grid_sizeY/2,:,t_end-1));
end

