%% computes stress tensor and pressure at a point (x,y).
% Geometry is the unit circle and the density function is (sin t, 1)

close all

x = 10;
y = 1;

% exact pressure and stress tensor computed using Mathematica
p_exact = @(x,y) 2*x*y/(x^2+y^2)^2;
ux_exact = @(x,y) x*y*(3*x^2-2*x^4-3*y^2+2*y^4)/(x^2+y^2)^4;
uy_exact = @(x,y) -(-4*x^6 + 3*y^4 + 6*x^2*y^2*(2*y^2-3)+x^4*(3+8*y^2))/(4*(x^2+y^2)^4);
vx_exact = @(x,y) (-3*x^4+6*x^2*(3-2*x^2)*y^2 - (3+8*x^2)*y^4 + 4*y^6)/(4*(x^2+y^2)^4);
vy_exact = @(x,y) x*(x-y)*y*(x+y)*(-3+2*x^2+2*y^2)/(x^2+y^2)^4;

grad_exact = @(x,y) [ux_exact(x,y), uy_exact(x,y); vx_exact(x,y), vy_exact(x,y)];
stress_v_exact = @(x,y) grad_exact(x,y) + grad_exact(x,y)';
stress_exact = @(x,y) stress_v_exact(x,y) - p_exact(x,y)*eye(2);

prams.nv = 1; % number of bodies
prams.T = 10; % time horizon
prams.m = 1; % number of time steps
prams.lengths = 1*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 2;
prams.Nbd = 0; %points on solid wall
prams.nbd = 0; %number of walls

options.farField = 'extensional';
options.saveData = true;
options.fileBase = 'stress_convergence';
options.append = false;
options.inear = true;
options.usePreco = false;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 1;
options.confined = false;

[options,prams] = initRigid2D(options,prams);

xc = [0;0];
tau = 0;

% compute stress and pressure using various grid sizes
N = [8,16,32,64,128];
stress1 = zeros(2,length(N));
stress2 = zeros(2,length(N));
p_approx = zeros(1,length(N));

for i = 1:length(N)
    prams.N = N(i); 
    Xfinal = rigid2D(options, prams, xc, tau);

    pp = post(['../output/data/',options.fileBase, '.mat']);

    % set density function
    theta = (0:prams.N-1)'*2*pi/prams.N;
    pp.etaF = [sin(theta);ones(length(theta),1)];
    
    p_approx(i) = pp.evaluatePressure(1, [x;y], 'fibers');
    [stress1(:,i), stress2(:,i)] = pp.evaluateStress(1, [x;y], 'fibers');
end

stress_approx = [stress1(:,end),stress2(:,end)];

% find ratio of computed quantities to exact quantities
fsv = (stress_approx + p_approx(end)*eye(2))./stress_v_exact(x,y);
fs = stress_approx./stress_exact(x,y);
fp = p_approx(end)/p_exact(x,y);