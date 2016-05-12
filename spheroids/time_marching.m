close all;

N = 128;
theta = (0:N-1)'*2*pi/N;

a = 2;
b = 0.5;
tau1(1) = pi/4;
tau2(1) = pi/2;

c1 = [2;-0.2];
c2 = [-1;0.2];

particle1 = calculate_ellipsoid(a, b, c1, tau1(1), N);
particle2 = calculate_ellipsoid(a, b, c2, tau2(1), N);

figure();
fill(particle1.x,particle1.y,'k');
hold on;
fill(particle2.x,particle2.y,'k');
xlim([-20,20]);
ylim([-10,10]);
axis equal
drawnow;

t = 0;
tf = 10;
dt = 0.1;

i = 2;
 
while (t < tf)
    
    t = t+dt;
    A = stokes_DLP_matrix_2_particles_power_force_free(particle1, particle2);
    
    
    u_inf = -[particle1.y; zeros(N,1); particle2.y; zeros(N,1); zeros(6,1)];
    
    eta = A\u_inf;
    
    U1 = eta(4*N+1:4*N+2);
    U2 = eta(4*N+3:4*N+4);
    omega1 = eta(4*N+5);
    omega2 = eta(4*N+6);
    
    tau1(i) = tau1(i-1) + omega1*dt;
    tau2(i) = tau2(i-1) + omega2*dt;
    c1 = c1 + U1*dt;
    c2 = c2 + U2*dt;
    
    particle1 = calculate_ellipsoid(a, b, c1, tau1(i), N);
    particle2 = calculate_ellipsoid(a, b, c2, tau2(i), N);

    hold off
    fill(particle1.x,particle1.y,'k');
    hold on;
    fill(particle2.x,particle2.y,'k');
    xlim([-20,20]);
    ylim([-10,10]);
    axis equal
    drawnow;
    
    i=i+1;
end

% xOmega = linspace(-6,6);
% yOmega = linspace(-6,6);
% [X,Y] = meshgrid(xOmega, yOmega, 50);
% 
% u = zeros([size(X),2]);
% for i = 1:size(X,1)
%     for j = 1:size(X,2)
%         
%         if (X(i,j)^2 + Y(i,j)^2 > 1.1 && (X(i,j)-2)^2 + (Y(i,j)-2)^2 >1.1)
%             u(i,j,:) =  0*u_inf_handle(X(i,j),Y(i,j)) + evaluate_stokes_DLP_2_particles(X(i,j), Y(i,j), ...
%                         particle1, particle2, eta_vector1, eta_vector2);
%         else
%             u(i,j,:) = nan;
%         end
%     end
% end
% 
% close all
% figure()
% contourf(X,Y,u(:,:,1))
% hold
% fill(particle1.x,particle1.y,'k');
% fill(particle2.x,particle2.y,'k');
% title('u', 'fontsize', 15);
% axis equal
% colorbar
% 
% figure
% contourf(X,Y,u(:,:,2))
% hold 
% fill(particle1.x,particle1.y,'k');
% fill(particle2.x,particle2.y,'k');
% title('v', 'fontsize', 15);
% axis equal
% colorbar
% 
% skip = 2;
% figure
% quiver(X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),u(1:skip:end,1:skip:end,1),u(1:skip:end,1:skip:end,2),2)
% hold
% fill(particle1.x,particle1.y,'k');
% fill(particle2.x,particle2.y,'k');
% axis equal