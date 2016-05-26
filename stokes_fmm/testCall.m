N = 2^14;
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

tic
mu = deny - 1i*denx;
zn = nx + 1i*ny;
dip1 = 0.25/pi*mu.*zn;
dip2 = 0.25/pi*(mu.*conj(zn) - conj(mu).*zn);
%dip3 = 0.25/pi*(-2)*1i*(denx.*nx + deny.*ny);
vel = stokesDLPfmm(dip1,dip2,x,y);
v1 = real(vel);
u1 = -imag(vel);
toc

u2 = zeros(N,1);
v2 = zeros(N,1);
tic
for k = 1:N
  ind = [(1:k-1) (k+1:N)];
  rx = x(k) - x(ind);
  ry = y(k) - y(ind);
  rho2 = rx.^2 + ry.^2;
  rdotn = rx.*nx(ind) + ry.*ny(ind);
  rdotden = rx.*denx(ind) + ry.*deny(ind);
  kernel = rdotn.*rdotden./(rho2.^2);
  u2(k) = sum(kernel.*rx);
  v2(k) = sum(kernel.*ry);
end
u2 = u2/pi;
v2 = v2/pi;
toc


norm(u1 - u2)/norm(u2)
norm(v1 - v2)/norm(v2)




