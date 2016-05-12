function particle = calculate_ellipsoid(a, b, xc, tau, N)
% draw_bodies(x, tau, a L) draws a set of rigid bodies with centres x,
% orientations tau, minor axis a and major axis L
%
% xc and tau are given as matrices of size Nx2, where N is the number of 
% particles. L and a are all vectors of length N.
%
% Algorithm taken from http://goo.gl/qCmkEj


theta = (0:N-1)'*2*pi/N;
x_ellipse = a*cos(theta) + xc(1);
y_ellipse = b*sin(theta) + xc(2);

particle.x = x_ellipse*cos(tau) - y_ellipse*sin(tau);
particle.y = x_ellipse*sin(tau) + y_ellipse*cos(tau);

xp = -a*sin(theta)*cos(tau) - b*cos(theta)*sin(tau);
yp = -a*sin(theta)*sin(tau) + b*cos(theta)*cos(tau);
xpp = -a*cos(theta)*cos(tau) + b*sin(theta)*sin(tau);
ypp = -a*cos(theta)*sin(tau) - b*sin(theta)*cos(tau);


particle.jac = sqrt(xp.^2 + yp.^2);
particle.tau_x = xp./particle.jac;
particle.tau_y = yp./particle.jac;

particle.n_x = particle.tau_y;
particle.n_y = -particle.tau_x;

cur_n = xp.*ypp - yp.*xpp;
cur_d = (xp.^2 + yp.^2).^(3/2);
particle.cur = -cur_n./cur_d;
particle.c = xc';
particle.N = N;