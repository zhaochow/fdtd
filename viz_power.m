function viz_power(evo_time)
% Ez and power along X axis
[grid_sizeX, grid_sizeY, t_end] = size(evo_time);
t_final = t_end -1;
Xaxis = 1:grid_sizeX;
Ez = zeros(1,grid_sizeX);
power = zeros(1,grid_sizeX);
for i = 1:grid_sizeX
    Ez(i) = rms(squeeze(evo_time(grid_sizeY/2,i,t_final-1/(freq*deltaT):t_final))*sqrt(2));
    power(i) = Ez(i)^2/(2*120*pi);
end
figure
plot(Xaxis,Ez);
title("Electric field")
figure
semilogy(Xaxis,power);
title("Power")
end

