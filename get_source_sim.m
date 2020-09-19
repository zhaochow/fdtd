function sim = get_source_sim()
    c = 299792458;                                                          % Speed of light
    grid_sizeX = 600;
    grid_sizeY = 500;
    freq = 500e12;
    lambda = c/freq;
    deltaX = lambda/20;                     %Half-boosted precision
    deltaT = deltaX/(c*sqrt(2));
    
    sim = get_base_sim(grid_sizeX, grid_sizeY, deltaX, deltaT);
    
    
    sim.source_idx = [grid_sizeY/2,grid_sizeX/2];
    sim.source_coeff = [2*pi*freq*deltaT];
    sim.source_phase = [0];
    
end

