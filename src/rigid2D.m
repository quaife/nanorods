function Xfinal = rigid2D(X,options,prams)

om = monitor(options,prams);
tt = tstep(options,prams);

oc = curve;
time = 0;
while time < prams.T
  time = time + tt.dt;
  geom = capsules(X);

  [density,Up,wp,iter,flag] = tt.timeStep(geom);

  [x,y] = oc.getXY(X);
  [cx,cy] = oc.getXY(geom.center);
  for k = 1:prams.nv
    x(:,k) = x(:,k) + tt.dt*(Up(1,k)*ones(prams.N,1) - ...
        wp(k)*(y(:,k) - cy(k)*ones(prams.N,1)));
    y(:,k) = y(:,k) + tt.dt*(Up(2,k)*ones(prams.N,1) + ...
        wp(k)*(x(:,k) - cx(k)*ones(prams.N,1)));
  end
  X = oc.setXY(x,y);

  om.outputInfo(X);
%  clf
%  plot(x,y,'r');
%  axis equal
%  axis(2*[-5 5 -5 5])
%  pause(1e-2)

end


Xfinal = X;


