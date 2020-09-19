function sim = get_base_sim(grid_sizeX, grid_sizeY, deltaX, deltaT)
    eps_vac = 8.85418782e-12;
    mu_vac = 4*pi*1e-7;
    c_mu = deltaT/(mu_vac*deltaX);
    c_eps = deltaT/(eps_vac*deltaX);
    
    sim = struct;
    
    sim.grid_sizeX = grid_sizeX;
    sim.grid_sizeY = grid_sizeY;
    sim.deltaX = deltaX;
    sim.deltaT = deltaT;
    sim.grid_unitY = grid_sizeY/10;
    sim.grid_unitX = grid_sizeX/10;
    
    sim.fields = zeros(grid_sizeY,grid_sizeX,3); % Ez(a,b,c),Hx(a,b+1/2,c),Hy(a+1/2,b,c)
    sim.param = zeros(grid_sizeY,grid_sizeX,3); % deltaT/(eps*deltaX),deltaT/(mu*deltaX)
    EPS = 1;
    MU = 2;
    SIG = 3;
    sim.param(:,:,EPS) = c_eps;
    sim.param(:,:,MU) = c_mu;
    sim.param(:,:,SIG) = 0;
end

