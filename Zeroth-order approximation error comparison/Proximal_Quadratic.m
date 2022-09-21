function Proximal_Quadratic

rng(60607)

n = 4;
iters = 30;

Alpha = 0.1; a = 0.1; % Factor of unknown nonliear function

Beta = 0.5;  % Factor of Line search

P = rand_psdmat(n);
Q = rand_condmat(n,10);

hat_J = @(x) (1/2*x'*P*x);

dhat_J = @(x) (P*x);

delta_J = @(x) (Alpha*1/2*(x-a)'*Q*(x-a));


% [x_equi,~] = Exact_Solution(hat_J,delta_J,n);

x_star =  Alpha * ( ( P + Alpha * Q ) \ Q ) * ones(n,1) * a;


J = @(x) ( hat_J(x) + delta_J(x) );


x0 = Alpha*a /( 1 + Alpha);
x_list = zeros(n,iters);
x_list(:,1) = x0;

total_count = 0;
x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;


for i = 1:iters-1
    
    last_count = total_count;
    x = x_list(:,i);

    delta = 1e-3;
    
    dJ = get_grad(J, x, delta);
    total_count = total_count + 2 * n + 1;
    
    grad  = dJ - dhat_J(x);
    
    t = bt_line_search_Proximal(hat_J, x, grad, Beta, P, n);
%     total_count = total_count + count;
    
    u  = x - t * grad;
    
    x_list(:,i+1) = ( P + 1/t * eye(n) ) \ ( 1/t * u );
    
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

% plot( 1:iters, norms(x_diff, [], 1) )







end