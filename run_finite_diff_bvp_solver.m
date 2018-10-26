%runs finite_diff_bvp_solver

errors=[];
x=linspace(0,1,1000); % true solution grid
u = @(x) x.^3 + cos(pi*x);

hold on
plot(x,u(x),'Linewidth',2);
n=10;
for i=1:5;
    h=1/n;
    [approx_sol,error]=finite_diff_bvp_solver(0,1,n);
    errors(i) = error;
    %x=[h:h:(n-1)*h]';
    x=[0:h:(n-1)*h]';
    plot(x,approx_sol,'Linewidth',2);
    n=2*n;
end
legend('true','10','20','40','80','160');
hold off
errors'
ratios = errors(1:4)./errors(2:5);
ratios'

    
    
    
