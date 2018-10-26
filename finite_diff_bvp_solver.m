function [approx_sol,error] = finite_diff_bvp_solver(a,b,n)

%set up function
r = @(x) exp(x);
f = @(x) ((-6.*x)+((pi.^2).*cos(pi.*x))+(exp(x).*(x.^3+cos(pi*x))));

%set up boundary conditions
alpha = 1;
beta = 0;
gamma = 0;

%set up the grid
h=(b-a)/n;
%x=[h:h:(n-1)*h]';
x=[0:h:(n-1)*h]';

%assemble the matrix A (store as sparse)
%A = sparse(diag(2 + h^2*r(x),0) + diag(-ones(n-2,1),1) + diag(-ones(n-2,1),-1));
A = sparse(diag(2 + h^2*r(x),0) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1));
A(1,2) = -2;
A = (1/h^2)*A;


%assemble right-hand side vector
b= f(x);

%add boundary conditions
% b(1) = b(1) + alpha/h^2; %problem 1
% b(n-1) = b(n-1) + beta/h^2; %problem 1
b(1) = b(1)-(2*gamma)/h;
b(n) = b(n) + beta/h^2;

%compute approximate solution
approx_sol = A\b;

%compute error
true_solution = x.^3 + cos(pi*x);
error = max(abs(true_solution - approx_sol));

%plot approx solution and true solution

end
