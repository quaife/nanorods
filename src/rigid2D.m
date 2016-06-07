function Xfinal = rigid2D(options, prams, xc, tau)

om = monitor(options, prams);
tt = tstep(options,prams);

om.writeData(0, xc, tau, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));

time = 0;

while time < prams.T
    
    tic;
    
    time = time + tt.dt;
    geom = capsules(prams, xc, tau);
    
    [density,Up,wp,iter,flag, res] = tt.timeStep(geom);
    
    xc = xc + tt.dt*Up; %update centres
    tau = tau + tt.dt*wp; %update angles
    X = geom.getXY();

    om.writeMessage(....
        ['Finished t=', num2str(time), ' in ' num2str(iter) ' iterations after ', num2str(toc), ' seconds (residual ', num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING GMRES flag = ', num2str(flag)]); 
    end
    
    om.writeData(time, xc, tau, Up(1,:), Up(2,:), wp);
    om.writeDensity(time, density);   
    
end

Xfinal = X;


