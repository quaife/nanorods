function A = stokes_DLP_matrix(x, y, nx, ny, tau_x, tau_y, cur, jac, c)
%% stokes_DLP_matrix Create matrix corresponding to the double layer
% potential for unbounded Stokes flow around a single object
% A = stokes_DLP_matrix(x, y, nx, ny, tau_x, tau_y, cur, jac, c) creates a
% matrix for flow around an object with boundary given by [x,y], normal
% vector [nx, ny], tangential vector [tau_x, tau_y], curvature cur,
% arclength jac and centre c. 


N = length(x);
A = zeros(2*N + 3, 2*N + 3);

% make matrix for double layer potential
for i = 1:N
   for j = 1:N 
       
       n_outer = [nx(i), ny(i)]'*[nx(j), ny(j)];
       
       if (i ~= j)
           r = [x(i), y(i)] - [x(j), y(j)];
           
           r_outer = r'*r;
           rdotn = r(1)*nx(j) + r(2)*ny(j);
           rho = norm(r);          
           
           A(i,j) = (1/pi)*r_outer(1,1)*rdotn*jac(j)/rho^4;% + jac(j)*n_outer(1,1);
           A(i,j+N) = (1/pi)*r_outer(1,2)*rdotn*jac(j)/rho^4;% + jac(j)*n_outer(1,2);
           
           A(i+N,j) = (1/pi)*r_outer(2,1)*rdotn*jac(j)/rho^4;% + jac(j)*n_outer(2,1);
           A(i+N,j+N) = (1/pi)*r_outer(2,2)*rdotn*jac(j)/rho^4;% + jac(j)*n_outer(2,2);
       else
           
           tau_outer = [tau_x(j), tau_y(j)]'*[tau_x(j), tau_y(j)];           

           
           A(i,j) = -(1/pi)*tau_outer(1,1)*cur(j)*jac(j)/2;% + jac(j)*n_outer(1,1);
           A(i,j+N) = -(1/pi)*tau_outer(1,2)*cur(j)*jac(j)/2;% + jac(j)*n_outer(1,2);
           
           A(i+N,j) = -(1/pi)*tau_outer(2,1)*cur(j)*jac(j)/2;% + jac(j)*n_outer(2,1);
           A(i+N,j+N) = -(1/pi)*tau_outer(2,2)*cur(j)*jac(j)/2;% + jac(j)*n_outer(2,2);
           
       end
       
   end
end


A = 2*pi/N*A - 0.5*eye(2*N+3,2*N+3);

% add Stokeslet and Rotlet terms
for i = 1:N
    r = [x(i), y(i)] - c;
    rho = norm(r);
    
    r_outer = r'*r;
    
    A(i,2*N+1) = 1/(4*pi)*(-log(rho) + r_outer(1,1)/rho^2);
    A(i+N,2*N+1) =  1/(4*pi)*(r_outer(2,1)/rho^2);
    
    A(i,2*N+2) = 1/(4*pi)*(r_outer(1,2)/rho^2);
    A(i+N,2*N+2) =  1/(4*pi)*(-log(rho) + r_outer(2,2)/rho^2);
    
    r_perp = [r(2), -r(1)];
    
    A(i,2*N+3) = r_perp(1)/rho^2;
    A(i+N,2*N+3) = r_perp(2)/rho^2;
    
end

for j = 1:N
   A(2*N+1,j) = (2*pi/N)*jac(j)/(2*pi);
   A(2*N+2,j+N) = (2*pi/N)*jac(j)/(2*pi);
   
   A(2*N+1, 2*N+1) = -1;
   A(2*N+2, 2*N+2) = -1;
   
   x_perp = [y(j), - x(j)];
   
   A(2*N+3,j) = (2*pi/N)*x_perp(1)*jac(j)/(2*pi);
   A(2*N+3,j+N) = (2*pi/N)*x_perp(2)*jac(j)/(2*pi);
   A(2*N+3,2*N+3) = -1;
   
end

