function Model_Free_Quadratic

rng(60607)

n = 5;
iters = 55;

a = 0.1; b = 0.1; % Factor of unknown nonliear function

rho = 0.3; Beta = 0.5;  % Factor of Line search

P = rand_psdmat(n);
Q = rand_condmat(n,10);

hat_f = @(x) (1/2*x'*P*x);

% nabla_hat_f = @(x) (P*x);

g = @(x) (a*1/2*(x-b)'*Q*(x-b));


% [x_equi,~] = Exact_Solution(hat_J,delta_J,n);

x_star =  a * ( ( P + a * Q ) \ Q ) * ones(n,1) * b;

f = @(x) ( hat_f(x) + g(x) );


x0 = -1e-7;
x_list = zeros(n,iters);
x_list(:,1) = x0;

total_count = 0;
x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;

for i = 1:iters-1
    
    x = x_list(:,i);
%     delta = 1e-3;
    last_count = total_count;
%     dJ = get_grad(J, x, delta);
    dJ = P * x + a * Q * ( x - b );
    total_count = total_count + 2 * n + 1;
    
    [eta,count] = bt_line_search_MF(f , x, dJ, rho, Beta);
    total_count = total_count + count;
    
    x_list(:,i+1) = x - eta * dJ;
    x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
    x_list_count(:, total_count) = x_list(:,i+1);
    
    
    
    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
        
end

% x_diff = x_list - x_star;
% 
% semilogy( 1:iters, norms(x_diff, [], 1))

x_list_count = x_list_count(:,1:total_count);
x_diff_count = x_list_count - x_star;

semilogy( 1:total_count, norms(x_diff_count, [], 1) )




end