function total_count = GC_Quad(Gamma)
rng(60607)

n = 9;
% n  = 1;
iters = 100;

a = 0.1; b = 0.1; % Factor of unknown nonliear function

rho = 0.3; Beta = 0.5;  % Factor of Line search

% Gamma = 0.2; % Factor of Relative Error, \beta in Note 0905

P = rand_psdmat(n);
Q = rand_condmat(n,10);

% P = 1;
% 
% Q = 1;

L_hat_f = max(eigs(P));
L_f = max(eigs(P + a*Q));
% keyboard

hat_f = @(x) ( 1/2 * x' * P * x );


d_hat_f = @(x) ( P * x );

r = @(x) ( a * 1/2 * (x-b)' * Q * (x-b) );

d_r = @(x) ( a * Q * (x-b) );


f = @(x) ( hat_f(x) + r(x) );

d_f = @(x) ( d_hat_f(x) + d_r(x) );

x_star =  a * ( ( P + a * Q ) \ Q ) * ones(n,1) * b;


x0 = 0;

x_list = zeros(n,iters);
x_list(:,1) = x0;

% State = 0: Model-Based Method
% State = 1: Model-Free Method

delta_G0  = zeros(n,1);

delta_G_list = zeros(n,100);

delta_G_list(:,1) = delta_G0;

total_count = 1;

x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;

MF_count = 1;
MB_count = 1;
new_state = 0;
L_r = L_f - L_hat_f;
E_x = 0;
MF_Start_Index = [];
MF_End_Index = [];
MB_Start_Index = [];
MB_End_Index = [];

tol = 1e-3;
i = 0;
x = x0; 
% for i = 1:iters-1
while norm(x - x_star) > tol
    i = i + 1;
    x = x_list(:,i);
    delta_G = delta_G_list(:,MF_count);
    
    grad = d_hat_f(x);
 
    State = new_state;
    
    if State == 0
       %% Model-Based Regime
        [eta, LS_MB_fail, count] = bt_LS_Relative_Shuo_MB(f, x, grad + delta_G, rho, Beta, State);

        x_list(:,i+1) = x - eta * (grad + delta_G);
        x_plus = x_list(:,i+1);
       
        disp('Model Base')
        
        
        
%         if norm( nabla_g(x) - delta_G ) > ( 1 - Gamma ) / Gamma * norm( nabla_hat_f(x) + delta_G )
        if E_x > ( 1 - Gamma ) / Gamma * norm( d_hat_f(x) + delta_G )
            large_rel_error = 1;
            disp('Large Relative Error')
        else
            large_rel_error = 0;
        end
        
        E_x_plus = ( L_f - L_hat_f ) * norm(x_plus - x) + E_x;


%         count = count + 2 * n + 1;
        if LS_MB_fail == 0 && large_rel_error == 0
            new_state = 0;
            E_x = E_x_plus;
        else
            new_state = 1;
%             count = 1;
        end
        
        last_count = total_count;
        total_count = total_count + count;         
        
        MB_Start_Index(MB_count) = last_count ;
        MB_End_Index(MB_count) = total_count;
        MB_count = MB_count + 1;
        
        
        x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
        x_list_count(:, total_count) = x_list(:,i+1);
        
        

    else
       %% Model-Free Regime
       
       
       disp('Model Free')

       % Compare with analytical solution
       dJ_real = P * x + a * Q * ( x - b );
       
       
       last_count = total_count;
       total_count = total_count + 2 * n + 1;
       
       [eta,count] = bt_line_search_MF(f, x, dJ_real, rho, Beta);
%        eta = 0.01; 
%        err = norm( dJ - dJ_real )
%        count = 0;
       
       total_count = total_count + count;
%        keyboard
       MF_Start_Index(MF_count) = last_count ;
       MF_End_Index(MF_count) = total_count;
       
       x_list(:,i+1) = x - eta * dJ_real;

       x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
       x_list_count(:, total_count) = x_list(:,i+1);
       delta_G_list(:,MF_count + 1) = dJ_real - grad;
       MF_count = MF_count + 1;
       new_state = 0;
       
       E_x = ( L_r ) * norm(x_list(:,i+1) - x);
       
       
    end
    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
        
end

% x_diff = x_list - x_star;


% x_list_count = x_list_count(:,1:total_count);
% x_diff_count = x_list_count - x_star;
% 
% 
% for i = 1: MF_count - 1
%     MF_steps = MF_Start_Index(i)  : MF_End_Index(i);
% %     keyboard
%     MF_regime = x_diff_count(:, MF_steps);
% %     MF_diff = reshape(MF_regime,[n*p,size(MF_regime,3)]);
% %     keyboard
%     Model_free_progeress = norms(MF_regime,[], 1);
%     
%     semilogy( MF_steps, Model_free_progeress , 'b','LineWidth',1.5)
%     hold on
% end
% 
% for i = 1: MB_count - 1
%     MB_steps = MB_Start_Index(i): MB_End_Index(i);
%     MB_regime = x_diff_count(:, MB_steps);
%     
% %     keyboard
%     Model_based_progeress = norms(MB_regime,[], 1);
%     
%     semilogy( MB_steps, Model_based_progeress , 'r','LineWidth',1.5)
%     hold on
% end
% % title('GC: Quadratic Function with $\gamma = 0.9$','Interpreter','latex','FontSize',16)
% xlabel('Number of Function Evaluations','FontSize',16)
% ylabel('$\|x - x^{*}\|$','Interpreter','latex','FontSize',16)
% 
% 
%     if Gamma == 0.9
%         print('GC_Quad_large_gamma','-depsc2')
%     end
%     
%     if Gamma == 0.6
%         print('GC_Quad','-depsc2')
%     end


end








