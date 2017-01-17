close all


prams.N = 64; % points per body
prams.Nbd = 0; %points on solid wall

g = -5;%shear rate

prams.nv = 1; % number of bodies
prams.nbd = 0; %number of walls
prams.lengths = 1;
prams.widths = 0.25;
prams.T = pi*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g)/2;
prams.m = 50; % number of time steps
prams.order = 2;

options.farField = 'shear';
options.saveData = true;
options.fileBase = 'jeffery';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.confined = false;
options.dt_min = 1e-5;
options.rk_max_up = 1.5;
options.rk_max_down = 0.5;
options.rk_safety = 0.9;

[options,prams] = initRigid2D(options,prams);
xWalls = [];

xc = [0;0];
tau  = pi/2;
rk_tol = [1e-3,1e-4,1e-5,1e-6,1e-7];
adaptive_orders = [4, 5];

errors_adaptive = zeros(length(rk_tol),length(adaptive_orders));
cpu_times_adaptive = zeros(length(rk_tol),length(adaptive_orders));
rejects_adaptive = zeros(length(rk_tol),length(adaptive_orders));

% define exact solutions
a = prams.widths/2;
b = prams.lengths/2;
tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);

times_approx = cell(length(adaptive_orders)*length(rk_tol));
tau_approx = cell(length(adaptive_orders)*length(rk_tol));

k = 1;
for i = 1:length(rk_tol)
    options.rk_tol = rk_tol(i);
    
    for j = 1:length(adaptive_orders)
    
        disp(['tol = ', num2str(rk_tol(i)), ', order = ',...
                                    num2str(adaptive_orders(j))]);
        
        options.tstep_order = adaptive_orders(j);
        tic;
        [~, rejects_adaptive(i,j)] = rigid2DAdaptive(options, prams, xc, tau, xWalls);
        cpu_times_adaptive(i,j) = toc;
        
        pp = post(['../output/data/',options.fileBase, '.mat']);
        errors_adaptive(i,j) = abs(tau_exact(pp.times(end-1)) + pi/2 - wrapTo2Pi(pp.tau(end-1)));
        
        times_approx{k} = pp.times;
        tau_approx{k} = pp.tau;
        
        k = k+1;
    end   
end

loglog(errors_adaptive, cpu_times_adaptive, '-o')
legend('RK23', 'RK45', 'RK54');
xlabel('global error');
ylabel('CPU time');