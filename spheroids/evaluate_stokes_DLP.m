function u = evaluate_stokes_DLP(xOmega, xBoundary, eta, n, jac)


u = zeros(2,1);
N = length(xBoundary);

for i = 1:N
    
    r = xOmega - xBoundary(i,:);
    rho = norm(r);
    r_outer = r'*r;
    
    u = u + (r_outer/rho^4)*eta(i,:)'*jac(i)*r*n(i,:)';
    
end

u = 2*u/N;

