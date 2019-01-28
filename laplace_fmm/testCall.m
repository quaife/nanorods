%addpath ../src
%Test the FMM for stokes SLP.
N = 2^12;
theta = (0:N-1)'*2*pi/N;

q = rand(N,1)/N;
% charge

xs = rand(N,1);
ys = rand(N,1);
%source/target locations
dir1 = rand(N,1);
dir2 = rand(N,1);


tfmm = tic;
pot = laplaceDLPfmm(q,xs,ys,dir1,dir2);
tfmm = toc(tfmm);

potExact = zeros(N,1);
tdirect = tic;
for i = 1:N
  for j = 1:N
    if i ~= j
      rdotn = (xs(i) - xs(j))*dir1(j) + (ys(i) - ys(j))*dir2(j);
      rho2 = (xs(i) - xs(j))^2 + (ys(i) - ys(j))^2;
      potExact(i) = potExact(i) - q(j)*rdotn/rho2;
    end
  end
end
tdirect = toc(tdirect);
potExact = potExact/2/pi;


fprintf('FMM    required %4.2e seconds\n',tfmm);
fprintf('Direct required %4.2e seconds\n',tdirect);
fprintf('Error is        %4.2e\n',norm(pot - potExact)/norm(pot));
