N = 128;
theta = (0:N-1)'*2*pi/N;

particle1.x = cos(theta);
particle1.y = sin(theta);
particle1.n_x = cos(theta);% points towards outside of circle
particle1.n_y = sin(theta);
particle1.tau_x = -sin(theta);
particle1.tau_y = cos(theta);
particle1.cur = -ones(N,1);
particle1.jac = ones(N,1);
particle1.c = [0,0];
particle1.N = N;

% a = 1;
% b = 0.5;
% particle1.x = a*cos(theta);
% particle1.y = b*sin(theta);
% 
% particle1.jac = sqrt(a^2*sin(theta).^2 + b^2*cos(theta).^2);
% particle1.n_x = b*cos(theta)./particle1.jac;
% particle1.n_y = a*sin(theta)./particle1.jac;
% cur_n = a*b;
% cur_d = (a^2*sin(theta).^2 + b^2*cos(theta).^2).^(3/2);
% particle1.cur = -cur_n./cur_d;
% particle1.c = [0,0];
% particle1.N = N;
% particle1.tau_x = -a*sin(theta)./particle1.jac;
% particle1.tau_y = b*cos(theta)./particle1.jac;

particle2.x = cos(theta)+4;
particle2.y = sin(theta)+4;
particle2.n_x = cos(theta);% points towards outside of circle
particle2.n_y = sin(theta);
particle2.tau_x = -sin(theta);
particle2.tau_y = cos(theta);
particle2.cur = -ones(N,1);
particle2.jac = ones(N,1);
particle2.c = [4,4];
particle2.N = N;

A = stokes_DLP_matrix_2_particles_power_force_free(particle1, particle2);


u_inf = -[particle1.y; zeros(N,1); particle2.y; zeros(N,1); zeros(8,1)];
u_inf_handle = @(x,y) [y;0];

U = 0.01;
V = 0;
u_inf = -[U*ones(N,1);V*ones(N,1);U*ones(N,1);V*ones(N,1);zeros(6,1)];
u_inf_handle = @(x,y) [U;V]; 

% xc = 0;
% yc = 0;
% u_inf = -[-(particle1.y-yc)./((particle1.x-xc).^2+(particle1.y-yc).^2);(particle1.x-xc)./((particle1.x-xc).^2+(particle1.y-yc).^2);...
%     -(particle2.y-yc)./((particle2.x-xc).^2+(particle2.y-yc).^2);(particle2.x-xc)./((particle2.x-xc).^2+(particle2.y-yc).^2);zeros(6,1)];
% u_inf_handle = @(x,y) [-(y-yc)./((x-xc).^2+(y-yc).^2);(x-xc)./((x-xc).^2+(y-yc).^2)]; 
u_inf = -[particle1.x; zeros(3*N+6,1)];
u_inf_handle = @(x,y) [U;V]; 

eta = A\u_inf;

eta_vector1 = [eta(1:N), eta(N+1:2*N)];
eta_vector2 = [eta(2*N+1:3*N), eta(3*N+1:4*N)];
U1 = eta(4*N+1:4*N+2);
U2 = eta(4*N+3:4*N+4);
omega1 = eta(4*N+5);
omega2 = eta(4*N+6);
beta = eta(4*N+7:end);

lambda1 = zeros(2,1);
lambda2 = zeros(2,1);
xi1 = 0;
xi2 = 0;


%% plot results
% 
% r = linspace(1.1,50);
% omega = linspace(0,2*pi,33);
% [r,omega] = meshgrid(r,omega);
% X = r.*cos(omega);
% Y = r.*sin(omega);

xOmega = linspace(-6,6);
yOmega = linspace(-6,6);
[X,Y] = meshgrid(xOmega, yOmega, 50);

u = zeros([size(X),2]);
for i = 1:size(X,1)
    for j = 1:size(X,2)
        
        if (X(i,j)^2 + Y(i,j)^2 > 1.1 && (X(i,j)-4)^2 + (Y(i,j)-4)^2 >1.1)
            u(i,j,:) =  evaluate_stokes_DLP_2_particles(X(i,j), Y(i,j), ...
                        particle1, particle2, eta_vector1, eta_vector2, lambda1, xi1, lambda2, xi2, beta);
        else
            u(i,j,:) = nan;
        end
    end
end

close all
figure()
contourf(X,Y,u(:,:,1))
hold
fill(particle1.x,particle1.y,'k');
fill(particle2.x,particle2.y,'k');
title('u', 'fontsize', 15);
axis equal
colorbar

figure
contourf(X,Y,u(:,:,2))
hold 
fill(particle1.x,particle1.y,'k');
fill(particle2.x,particle2.y,'k');
title('v', 'fontsize', 15);
axis equal
colorbar

skip = 2;
figure
quiver(X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),u(1:skip:end,1:skip:end,1),u(1:skip:end,1:skip:end,2),2)
hold
fill(particle1.x,particle1.y,'k');
fill(particle2.x,particle2.y,'k');
axis equal



