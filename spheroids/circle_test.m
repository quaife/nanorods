N = 128;
theta = (0:N-1)'*2*pi/N;

x = cos(theta);
y = sin(theta);
nx = -cos(theta);% points towards inside of circle
ny = -sin(theta);

tau_x = sin(theta);
tau_y = -cos(theta);

cur = ones(N,1);

jac = ones(length(nx),1);
c = [0,0];

A = stokes_DLP_matrix(x, y, nx, ny, tau_x, tau_y, cur, jac, c);

U = 0.01;
u_inf = [U*ones(N,1);zeros(N+3,1)];
eta = A\u_inf;

%eta_vector = [eta(1:2:end-1), eta(2:2:end)];
eta_vector = [eta(1:N), eta(N+1:2*N)];

%% plot results

xOmega = linspace(-3,3);
[X,Y] = meshgrid(xOmega,xOmega, 200);

xBoundary = [x, y];
n = [nx,ny];

u = zeros([size(X),2]);
for i = 1:length(xOmega)
    for j = 1:length(xOmega)
        
        if (X(i,j)^2 + Y(i,j)^2 > 1.1)
            u(i,j,:) = [U;0] + evaluate_stokes_DLP([X(i,j),Y(i,j)], xBoundary, eta_vector, n, jac);
        else
            u(i,j,:) = nan;
        end
    end
end

close all
figure
contourf(X,Y,u(:,:,1))
axis equal
colorbar

figure
contourf(X,Y,u(:,:,2))
axis equal
colorbar

figure
quiver(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),u(1:4:end,1:4:end,1),u(1:4:end,1:4:end,2),2)
axis equal



