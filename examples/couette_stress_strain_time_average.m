rOut = 10;
rIn = 5;

strain = @(w) 2*w*rOut^2/(rOut^2-rIn^2);
omega = 1:5;

Nparticles = [5, 10, 20];
Nruns = 3;

runs = cell(length(Nparticles),1);

prams.N = 8; % points per body
prams.Nbd = 192; %points on solid wall

prams.nv = 2; % number of bodies
prams.nbd = 2; %number of walls
prams.T = 10*pi; % time horizon
prams.m = 500; % number of time steps
prams.lengths = 0.1*ones(1, prams.nv);
prams.widths = 0.1*ones(1,prams.nv);
prams.order = 2;
prams.tracker_fnc = @(t) [10,0;5*cos(t),5*sin(t)];
prams.gmresTol = 1e-8;

options.farField = 'couette';
options.saveData = true;
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = true;
options.rk_tol = 1e-2;
options.dt_min = 1e-5;
options.rk_max_up = 1.5;
options.rk_max_down = 0.25;
options.rk_safety = 0.9;

[options,prams] = initRigid2D(options,prams);


for i = 1:length(Nparticles)
    
    
    N = Nparticles(i);
    
    runs{i}.N = N;
    runs{i}.volFrac = zeros(Nruns,1);
    runs{i}.nv = zeros(Nruns,1);
    runs{i}.seeds = zeros(Nruns,1);
    runs{i}.omega = zeros(length(omega),1);
    runs{i}.stress = zeros(Nruns,prams.m,length(omega)); 
    runs{i}.strain = zeros(length(omega),1);
    
    for k = 1:Nruns
        c = clock;
        seed = c(end);        
        runs{i}.seeds(k) = seed;
        
        % add particles randomly
        prams.nv = 2;
        
        xc = [5, 5;-5, 5];
        tau = [0, 0];

        geom = capsules(prams, xc, tau);
        [xc, tau] = geom.fill_couette(5.26, 9.74, N, prams, seed);
        
        prams.nv = length(tau);        
        runs{i}.volFrac(k) = prams.nv*(prams.lengths(1)/2)^2/(10^2-5^2);
        runs{i}.nv = prams.nv;
        
        for j = 1:length(omega)
            
            options.couetteSpeed = omega(j);
            
            runs{i}.omega(j) = omega(j);
            runs{i}.strain(j) = strainRate(j);
            
            disp(['i = ', num2str(i), ' k = ', num2str(k), ' j = ', num2str(j)]);
            couette;
            
            for m = 1:prams.m
                runs{i}.stress(k,m,j) = norm(pp.stokeslets(:,:,m));
            end
            
            runs{i}.times = pp.times;
            save('stressDataTA', 'runs');
            %[~,stress_xy(k,i,j),~,~] = pp.computeWallStress(1);
        end
    end
end