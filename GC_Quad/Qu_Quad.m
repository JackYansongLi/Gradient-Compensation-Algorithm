function Qu_Quad
rng(60607)

n = 9;
iters = 100;

c2 = 0.1; c3 = 0.1; % Factor of unknown nonlinear function

Alpha = 0.3; Beta = 0.5;  % Factor of Line search

P = rand_psdmat(n);
Q = rand_condmat(n,10);


hat_f = @(x) ( 1/2 * x' * P * x );
r = @(x) ( c2 * 1/2 * (x-c3)' * Q * (x-c3) );
x_star =  c2 * ( ( P + c2 * Q ) \ Q ) * ones(n,1) * c3;
f = @(x) ( hat_f(x) + r(x) );


x0 = 0;

x_list = zeros(n,iters);
x_list(:,1) = x0;

% State = 0: Model-Based Method
% State = 1: Model-Free Method

total_count = 0;
count = 0;

x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;

MF_count = 1;


tol = 1e-3;
i = 0;
x = x0; 
% for i = 1:iters-1
while norm(x - x_star) > tol
    i = i + 1;
    x = x_list(:,i);
    last_count = total_count;
    total_count = total_count + count; 

    dJ_real = P * x + c2 * Q * ( x - c3 );



    total_count = total_count + 2 * n + 1;

    [eta,count] = bt_line_search_MF(f , x, dJ_real, Alpha, Beta);

    total_count = total_count + count;



    x_list(:,i+1) = x - eta * dJ_real;

    x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
    x_list_count(:, total_count) = x_list(:,i+1);


    MF_count = MF_count + 1;

    disp('Model Free')

    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
        
end


x_list_count = x_list_count(:,1:total_count);
x_diff_count = x_list_count - x_star;
figure
semilogy( 1:total_count, norms(x_diff_count, [], 1),'LineWidth',1.5)

% title('Qu et al.: Quadratic Function','FontSize',16)
xlabel('Number of Function Evaluations','FontSize',16)
ylabel('$\|x - x^{*}\|$','Interpreter','latex','FontSize',16)

print('Qu_Quad','-depsc2')

end