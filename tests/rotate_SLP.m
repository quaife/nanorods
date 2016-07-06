addpath ../src

N = 32;

theta = (0:N-1)'*2*pi/N;
x1 = cos(theta);
y1 = 2*sin(theta);

omega = 0.63;
x2 = cos(omega)*x1 - sin(omega)*y1;
y2 = sin(omega)*x1 + cos(omega)*y1;

geom1 = capsules([],[x1;y1]);
geom2 = capsules([],[x2;y2]);

om.profile = false;
op1 = poten(geom1,om);
op2 = poten(geom2,om);

D1 = op1.stokesDLmatrix(geom1);
D2 = op2.stokesDLmatrix(geom2);
% double layer potentials for the two geometries

etax = exp(sin(theta));
etay = sin(cos(exp(sin(theta))));

Deta1 = D1*[etax;etay];
Deta2 = D2*[etax;etay];


mux = cos(omega)*etax + sin(omega)*etay;
muy = -sin(omega)*etax + cos(omega)*etay;
Dmu = D1*[mux;muy];
Dmux = Dmu(1:N);
Dmuy = Dmu(N+1:end);
Deta3x = cos(omega)*Dmux - sin(omega)*Dmuy;
Deta3y = sin(omega)*Dmux + cos(omega)*Dmuy;
Deta3 = [Deta3x;Deta3y];

% norm of Deta2-Deta3 should be 0



