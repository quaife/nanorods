function A = stokes_DLP_matrix_2_particles_power_force_free(particle1, particle2)
%% stokes_DLP_matrix Create matrix corresponding to the double layer
% potential for unbounded Stokes flow around a single object
% A = stokes_DLP_matrix(x, y, nx, ny, tau_x, tau_y, cur, jac, c) creates a
% matrix for flow around an object with boundary given by [x,y], normal
% vector [nx, ny], tangential vector [tau_x, tau_y], curvature cur,
% arclength jac and centre c. 


N1 = particle1.N;
N2 = particle2.N;

N = N1+N2;
A = zeros(2*N + 6, 2*N + 6);

% make matrix for double layer potential

%% paricle 1
for i = 1:N1
    
    jac = particle1.jac;
    cur = particle1.cur;
    tau_x = particle1.tau_x;
    tau_y = particle1.tau_y;
    n_x = particle1.n_x;
    n_y = particle1.n_y;
    
    for j = 1:N1
        
        if (i ~= j)
            r = [particle1.x(i), particle1.y(i)] - [particle1.x(j), particle1.y(j)];
            
            r_outer = r'*r;
            rdotn = r(1)*n_x(j) + r(2)*n_y(j);
            rho = norm(r);
            
            A(i,j) = (1/pi)*r_outer(1,1)*rdotn*jac(j)/rho^4;
            A(i,j+N1) = (1/pi)*r_outer(1,2)*rdotn*jac(j)/rho^4;
            
            A(i+N1,j) = (1/pi)*r_outer(2,1)*rdotn*jac(j)/rho^4;
            A(i+N1,j+N1) = (1/pi)*r_outer(2,2)*rdotn*jac(j)/rho^4;
        else
            
            tau_outer = [tau_x(j), tau_y(j)]'*[tau_x(j), tau_y(j)];            
            
            A(i,j) = (1/pi)*tau_outer(1,1)*cur(j)*jac(j)/2;
            A(i,j+N1) = (1/pi)*tau_outer(1,2)*cur(j)*jac(j)/2;
            
            A(i+N1,j) = (1/pi)*tau_outer(2,1)*cur(j)*jac(j)/2;
            A(i+N1,j+N1) = (1/pi)*tau_outer(2,2)*cur(j)*jac(j)/2;
            
        end
        
    end
    
    jac = particle2.jac;
    n_x = particle2.n_x;
    n_y = particle2.n_y;
    for j = 1:N2
        
            r = [particle1.x(i), particle1.y(i)] - [particle2.x(j), particle2.y(j)];
            
            r_outer = r'*r;
            rdotn = r(1)*n_x(j) + r(2)*n_y(j);
            rho = norm(r);
            
            A(i,j+2*N1) = (1/pi)*r_outer(1,1)*rdotn*jac(j)/rho^4;
            A(i,j+2*N1+N2) = (1/pi)*r_outer(1,2)*rdotn*jac(j)/rho^4;
            
            A(i+N1,j+2*N1) = (1/pi)*r_outer(2,1)*rdotn*jac(j)/rho^4;
            A(i+N1,j+2*N1+N2) = (1/pi)*r_outer(2,2)*rdotn*jac(j)/rho^4;
    end
end

%% particle 2
for i = 1:N2
    
    jac = particle2.jac;
    cur = particle2.cur;
    tau_x = particle2.tau_x;
    tau_y = particle2.tau_y;
    n_x = particle2.n_x;
    n_y = particle2.n_y;
    
    for j = 1:N2
        
        if (i ~= j)
            r = [particle2.x(i), particle2.y(i)] - [particle2.x(j), particle2.y(j)];
            
            r_outer = r'*r;
            rdotn = r(1)*n_x(j) + r(2)*n_y(j);
            rho = norm(r);
            
            A(i+2*N1,j+2*N1) = (1/pi)*r_outer(1,1)*rdotn*jac(j)/rho^4;
            A(i+2*N1,j+2*N1+N2) = (1/pi)*r_outer(1,2)*rdotn*jac(j)/rho^4;
            
            A(i+2*N1+N2,j+2*N1) = (1/pi)*r_outer(2,1)*rdotn*jac(j)/rho^4;
            A(i+2*N1+N2,j+2*N1+N2) = (1/pi)*r_outer(2,2)*rdotn*jac(j)/rho^4;
        else
            
            tau_outer = [tau_x(j), tau_y(j)]'*[tau_x(j), tau_y(j)];            
            
            A(i+2*N1,j+2*N1) = (1/pi)*tau_outer(1,1)*cur(j)*jac(j)/2;
            A(i+2*N1,j+2*N1+N2) = (1/pi)*tau_outer(1,2)*cur(j)*jac(j)/2;
            
            A(i+2*N1+N2,j+2*N1) = (1/pi)*tau_outer(2,1)*cur(j)*jac(j)/2;
            A(i+2*N1+N2,j+2*N1+N2) = (1/pi)*tau_outer(2,2)*cur(j)*jac(j)/2;
            
        end
        
    end
    
    jac = particle1.jac;
    n_x = particle1.n_x;
    n_y = particle1.n_y;
    for j = 1:N2
        
            r = [particle2.x(i), particle2.y(i)] - [particle1.x(j), particle1.y(j)];
            
            r_outer = r'*r;
            rdotn = r(1)*n_x(j) + r(2)*n_y(j);
            rho = norm(r);
            
            A(i+2*N1,j) = (1/pi)*r_outer(1,1)*rdotn*jac(j)/rho^4;
            A(i+2*N1,j+N1) = (1/pi)*r_outer(1,2)*rdotn*jac(j)/rho^4;
            
            A(i+2*N1+N2,j) = (1/pi)*r_outer(2,1)*rdotn*jac(j)/rho^4;
            A(i+2*N1+N2,j+N1) = (1/pi)*r_outer(2,2)*rdotn*jac(j)/rho^4;
    end
end

A = -2*pi/N1*A - [0.5*eye(2*N,2*N), zeros(2*N,6); zeros(6,2*N+6)];

%% unknown velocity columns
for i = 1:N1
    A(i,2*N+1) = -1;
    A(i+N1,2*N+2) = -1;
    
    A(i+2*N1,2*N+3) = -1;
    A(i+2*N1+N2, 2*N+4) = -1;
    
    r = [particle1.x(i), particle1.y(i)] - particle1.c;
    A(i, 2*N+5) = -r(1);
    A(i+N1,2*N+5) = -r(2);
    
    r = [particle2.x(i), particle2.y(i)] - particle2.c;
    A(i+2*N1, 2*N+6) = -r(1);
    A(i+2*N1+N2,2*N+6) = -r(2); 
    
%     A(i,2*N+7) = 1;
%     A(i+N1, 2*N+8) = 1;
%     
%     A(i+2*N1,2*N+7) = 1;
%     A(i+2*N1+N2,2*N+8) = 1;
end

%% force free constraint
for j = 1:N1
    A(2*N+1,j) = 1*particle1.jac(j)/N1;
    A(2*N+2,j+N1) = 1*particle1.jac(j)/N1;
    
    A(2*N+3,j+2*N1) = 1*particle2.jac(j)/N1;
    A(2*N+4,j+2*N1+N2) = 1*particle2.jac(j)/N1;
end

%% torque free constraints
for j = 1:N1
   
   r = [particle1.x(j), particle1.y(j)] - particle1.c;
   A(2*N+5,j) =  -r(2)*particle1.jac(j)/N1;
   A(2*N+5,j+N1) = r(1)*particle1.jac(j)/N1;
   
   r = [particle2.x(j), particle2.y(j)] - particle2.c; 
   A(2*N+6, j+2*N1) = -r(2)*particle2.jac(j)/N1;
   A(2*N+6, j+2*N1+N2) = r(1)*particle2.jac(j)/N1;
end

%% beta constraints
for j = 1:N1
%    A(2*N+7,j) = 1*particle1.jac(j);
%    %A(2*N+7,j+2*N1) = 1*particle2.jac(j);
%    
%    %A(2*N+8,j + N1) = 1*particle1.jac(j);
%    A(2*N+8,j+2*N1+N2) = 1*particle2.jac(j);
%    
%    A(2*N+7,2*N+7) = -1;
%    A(2*N+8,2*N+8) = -1;
   
end