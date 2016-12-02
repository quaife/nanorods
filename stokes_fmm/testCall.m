close all

N_all = 2.^[4:14];
tests = 100;

time_fmm = zeros(length(N_all),1);
time_direct = zeros(length(N_all),1);

%N = 2^14;

for i = 1:length(N_all)
    N = N_all(i);
    disp(N);
    theta = (0:N-1)'*2*pi/N;
    
    x = rand(N,1);
    y = rand(N,1);
    nx = rand(N,1);
    ny = rand(N,1);
    denx = rand(N,1);
    deny = rand(N,1);
    
    %x = cos(theta);
    %y = sin(theta);
    %nx = cos(theta);
    %ny = sin(theta);
    %denx = cos(theta)*2*pi/N;
    %deny = cos(theta)*2*pi/N;
    
    for j = 1:tests
       
        mu = deny - 1i*denx;
        zn = nx + 1i*ny;
        dip1 = 0.25/pi*mu.*zn;
        dip2 = 0.25/pi*(mu.*conj(zn) - conj(mu).*zn);
        tic;
        %dip3 = 0.25/pi*(-2)*1i*(denx.*nx + deny.*ny);
        vel = stokesDLPfmm(dip1,dip2,x,y);
        v1 = real(vel);
        u1 = -imag(vel);
        time_fmm(i) = time_fmm(i) +  toc/tests;
        
%         u2 = zeros(N,1);
%         v2 = zeros(N,1);
%         tic
%         for k = 1:N
%             ind = [(1:k-1) (k+1:N)];
%             rx = x(k) - x(ind);
%             ry = y(k) - y(ind);
%             rho2 = rx.^2 + ry.^2;
%             rdotn = rx.*nx(ind) + ry.*ny(ind);
%             rdotden = rx.*denx(ind) + ry.*deny(ind);
%             kernel = rdotn.*rdotden./(rho2.^2);
%             u2(k) = sum(kernel.*rx);
%             v2(k) = sum(kernel.*ry);
%         end
%         u2 = u2/pi;
%         v2 = v2/pi;
        %time_direct(i) = time_direct(i) + toc/tests;
    end
    
end

loglog(N_all(1:length(time_direct)), time_direct, '-bo', 'linewidth',2);
hold on
loglog(N_all, time_fmm, '-ro', 'linewidth', 2);
xlabel('$N$', 'interpreter', 'latex');
ylabel('time (s)');
legend({'direct summation', 'FMM'}, 'location', 'NW');

% norm(u1 - u2)/norm(u2)
% norm(v1 - v2)/norm(v2)




