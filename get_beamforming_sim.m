function sim = get_beamforming_sim()
    c = 299792458;                                        % Speed of light
    grid_sizeX = 800;
    grid_sizeY = 800;
    freq = 2e9;
    lambda = c/freq;
    deltaX = lambda/10; 
    deltaT = deltaX/(c*sqrt(2));
    
    sim = get_base_sim(grid_sizeX, grid_sizeY, deltaX, deltaT);
    
    d = lambda/2;
    beta = 2*pi/lambda;
    rad_direction = 3*pi/4;
    d_pix = round(d/deltaX)
    phase_shift = beta*d*cos(rad_direction);
    sim.source_idx = [grid_sizeY/2,grid_sizeX/2;...
                      grid_sizeY/2,grid_sizeX/2-d_pix;...
                      grid_sizeY/2,grid_sizeX/2-2*d_pix;...
                      grid_sizeY/2,grid_sizeX/2-3*d_pix;...
                      ];
                  
    sim.source_coeff = [2*pi*freq*deltaT;...
                        2*pi*freq*deltaT;...
                        2*pi*freq*deltaT;...
                        2*pi*freq*deltaT;...
                        ];
                    
    sim.source_phase = [0;...
                        -phase_shift;...
                        -2*phase_shift;...
                        -3*phase_shift;...
                        ];
end
