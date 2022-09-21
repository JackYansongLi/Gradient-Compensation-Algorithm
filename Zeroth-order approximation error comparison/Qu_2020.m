function Qu_2020
rng(60607)

n = 4;
% n  = 1;
iters = 80;

c2 = 0.1; c3 = 0.1; % Factor of unknown nonlinear function

Alpha = 0.3; Beta = 0.5;  % Factor of Line search

P = rand_psdmat(n);
Q = rand_condmat(n,10);

% P = 1;
% 
% Q = 1;

hat_J = @(x) ( 1/2 * x' * P * x );

dhat_J = @(x) ( P * x );

delta_J = @(x) ( c2 * 1/2 * (x-c3)' * Q * (x-c3) );


% [x_equi,~] = Exact_Solution(hat_J,delta_J,n);

x_star =  c2 * ( ( P + c2 * Q ) \ Q ) * ones(n,1) * c3;





J = @(x) ( hat_J(x) + delta_J(x) );

% J_normal = @(x) ( hat_J(x) + a * 1/2 * (x-b)' * Q * x - a * 1/2 * x' * Q * b );


% x0 = randn(n,1);

% x0 = a*b /( 1 + a);

x0 = 0;

x_list = zeros(n,iters);
x_list(:,1) = x0;

% State = 0: Model-Based Method
% State = 1: Model-Free Method

% State = 0;
% delta_G0  = zeros(n,1);
% 
% delta_G_list = zeros(n,100);
% 
% delta_G_list(:,1) = delta_G0;

total_count = 0;
count = 0;

x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;

MF_count = 1;

for i = 1:iters-1
    x = x_list(:,i);
%     delta_G = delta_G_list(:,MF_count);
    
%     grad = dhat_J(x);
%     J_tilde = @(x) ( hat_J(x) + delta_G' * x );
    
%     [eta, new_state, count] = bt_line_search_MB(J_tilde , x, grad + delta_G, rho, Beta, State);
    last_count = total_count;
    total_count = total_count + count; 
%     State = new_state;
%     
%     if State == 0
%        %% Model-Based Method
%         x_list(:,i+1) = x - eta * (grad + delta_G);
% 
%         x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
%         x_list_count(:, total_count) = x_list(:,i+1);
%         
%         
%         disp('Model Base')
%     else
       %% Model-Free Method
       
       
%        delta = 1e-3;              
%        dJ = get_grad(J, x, delta);

       % Compare with analytical solution
       dJ_real = P * x + c2 * Q * ( x - c3 );
       
       

       total_count = total_count + 2 * n + 1;
       
       [eta,count] = bt_line_search_MF(J , x, dJ_real, Alpha, Beta);
%        eta = 0.01; 
%        err = norm( dJ - dJ_real )
       total_count = total_count + count;
       
       
       
       x_list(:,i+1) = x - eta * dJ_real;

       x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
       x_list_count(:, total_count) = x_list(:,i+1);
       
%        delta_G_list(:,MF_count + 1) = dJ_real - grad;
       MF_count = MF_count + 1;
       
       disp('Model Free')
%     end
    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
        
end

% x_diff = x_list - x_star;


x_list_count = x_list_count(:,1:total_count);
x_diff_count = x_list_count - x_star;
figure
semilogy( 1:total_count, norms(x_diff_count, [], 1) )

title('Finite difference')
xlabel('Number of Function Evaluations')
ylabel('$\|x - x^{*}\|$','Interpreter','latex')


% delta_G_list = delta_G_list(:,1:MF_count);
% 
% delta_G_star = a * Q * (x_star - b);
% 
% delta_G_diff  = delta_G_list - delta_G_star ;
% 
% WGCA = norms(delta_G_diff, [], 1)



% semilogy( 1:iters, norms(x_diff, [], 1) )
% 
% hold on

% plot( 1:iters, norms(x_diff, [], 1) )

end