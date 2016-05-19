function Xfinal = rigid2D(options,prams, varargin)

om = monitor(options,prams);
tt = tstep(options,prams);

om.writeData(0, varargin{1}, varargin{2}, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));

time = 0;
while time < prams.T
    
    tic;
    
    time = time + tt.dt;
    geom = capsules(prams, varargin{:});
    
    [density,Up,wp,iter,flag, res] = tt.timeStep(geom);
    
    varargin{1} = varargin{1} + tt.dt*Up; %update centres
    varargin{2} = varargin{2} + tt.dt*wp; %update angles
    X = geom.getXY();

    disp(['Finished t=', num2str(time), ' in ' num2str(iter) ' iterations after ', num2str(toc), ' seconds (residual ', num2str(res), ')']);
    
    om.writeData(time, varargin{1}, varargin{2}, Up(1,:), Up(2,:), wp);
    om.writeDensity(time, density);
end


Xfinal = X;


