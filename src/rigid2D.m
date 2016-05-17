function Xfinal = rigid2D(options,prams, varargin)

om = monitor(options,prams);
tt = tstep(options,prams);

om.writeData(0, varargin{1}, varargin{2});

time = 0;
while time < prams.T
    time = time + tt.dt;
    geom = capsules(prams, varargin{:});
    
    [density,Up,wp,iter,flag] = tt.timeStep(geom);
    
    varargin{1} = varargin{1} + tt.dt*Up; %update centres
    varargin{2} = varargin{2} + tt.dt*wp; %update angles
    X = geom.getXY();

    om.writeData(time, varargin{1}, varargin{2});    
end


Xfinal = X;


