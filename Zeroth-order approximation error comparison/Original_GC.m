function Original_GC
rng(60607)

n = 9;
% n  = 1;
iters = 170;

a = 0.1; b = 0.1; % Factor of unknown nonliear function

rho = 0.3; Beta = 0.5;  % Factor of Line search

P = rand_psdmat(n);
Q = rand_condmat(n,10);

% P = 1;
% 
% Q = 1;

hat_f = @(x) ( 1/2 * x' * P * x );

nabla_hat_f = @(x) ( P * x );

g = @(x) ( a * 1/2 * (x-b)' * Q * (x-b) );

nabla_g = @(x) ( a * Q * (x-b) );



% [x_equi,~] = Exact_Solution(hat_J,delta_J,n);

x_star =  a * ( ( P + a * Q ) \ Q ) * ones(n,1) * b;



f = @(x) ( hat_f(x) + g(x) );

% J_normal = @(x) ( hat_J(x) + a * 1/2 * (x-b)' * Q * x - a * 1/2 * x' * Q * b );


% x0 = randn(n,1);

% x0 = a*b /( 1 + a);

x0 = 0;

x_list = zeros(n,iters);
x_list(:,1) = x0;

% State = 0: Model-Based Method
% State = 1: Model-Free Method

State = 0;
delta_G0  = zeros(n,1);

delta_G_list = zeros(n,100);

delta_G_list(:,1) = delta_G0;

total_count = 0;

x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;

MF_count = 1;

% total_count = 1
new_state = 0;
for i = 1:iters-1
    x = x_list(:,i);
    delta_G = delta_G_list(:,MF_count);
    
    grad = nabla_hat_f(x);
    
    

    State = new_state;
    
    if State == 0
       %% Model-Based Method
        disp('Model Base')       
        [eta, new_state, count] = bt_line_search_MB(f , x, grad + delta_G, rho, Beta, State);
        last_count = total_count;
        total_count = total_count + count; 
        
        x_list(:,i+1) = x - eta * (grad + delta_G);

        x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
        x_list_count(:, total_count) = x_list(:,i+1);
        
        

    else
       %% Model-Free Method
%        keyboard;
       
       

       % Compare with analytical solution
       dJ_real = P * x + a * Q * ( x - b );
       
       

       total_count = total_count + 2 * n + 1;
       
       [eta,count] = bt_line_search_MF(f , x, dJ_real, rho, Beta);
%        eta = 0.01; 
%        err = norm( dJ - dJ_real )

       total_count = total_count + count;
       
       
       
       x_list(:,i+1) = x - eta * dJ_real;

       x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
       x_list_count(:, total_count) = x_list(:,i+1);
       
       delta_G_list(:,MF_count + 1) = dJ_real - grad;
       MF_count = MF_count + 1;
       new_state = 0;
       disp('Model Free')
    end
    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
        
end

% x_diff = x_list - x_star;


x_list_count = x_list_count(:,1:total_count);
x_diff_count = x_list_count - x_star;

figure;
semilogy( 1:total_count, norms(x_diff_count, [], 1) )

delta_G_list = delta_G_list(:,1:MF_count);

delta_G_star = a * Q * (x_star - b);

delta_G_diff  = delta_G_list - delta_G_star ;

norms(delta_G_diff, [], 1)


% semilogy( 1:iters, norms(x_diff, [], 1) )
% 
% hold on

% plot( 1:iters, norms(x_diff, [], 1) )

end








