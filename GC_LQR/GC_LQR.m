function GC_LQR

rng(60607)
iters = 100; 
Alpha = 0.3; Beta = 0.5; eta_min = 0.05;% Factor of Line Search
T = 50;                  % Factor of Function Evaluation, i.e., hat_C

r = 0.001;               % Gradient Estimation radius

Gamma = 0.6; L_r = 0.9;

n = 4; p = 3;
Q = 2*eye(n); R = eye(p);
A = randn(n,n);  B = randn(n,p); N = zeros(n,p);

[K_lin] = dlqr(A,B,Q,R,N);
mf_factor = 0.01;

f = @(x) [mf_factor*x(1)/(1-0.9*sin(x(1))); mf_factor*x(2)/(1-0.9*sin(x(2))); ...
          mf_factor*x(3)/(1-0.9*sin(x(3))); mf_factor*x(4)/(1-0.9*sin(x(4)))];
      

      
% State = 0: Model-Based Regime
% State = 1: Model-Free Regime

MF_count = 1;
MB_count = 1;
new_state = 0;
E_x = 0;
total_count = 1;
delta_G0  = zeros(p,n);
delta_G_list = zeros(p,n,100);
delta_G_list(:,:,1) = delta_G0;
K_list = zeros(p,n,iters);
K_list(:,:,1) = K_lin;
K_list_count = zeros(p,n,1000);
K_list_count(:,:,1) = K_lin;

x0 = randn(n,n);
mu = mean(x0,2);
mu2 = mu*mu';
sigma = cov(x0,1); 

MF_Start_Index = [];
MF_End_Index = [];
MB_Start_Index = [];
MB_End_Index = [];


K = K_lin;
load('K_star.mat')
i = 0;
tol = 1e-4;
% for i = 1:iters-1
while norm(K - K_star, 'fro') > tol
    i = i+1;
    K = K_list(:,:,i);
    
    delta_G = delta_G_list(:,:,MF_count);
    
    Pk = dlyap( (A-B*K)' , Q + K'*R*K );
    Ek = (R+B'*Pk*B)*K-B'*Pk*A;
    sigma_k = dlyap((A-B*K),sigma+mu2);
    grad = 2 * Ek * sigma_k;

    
    
    State = new_state;

    if State == 0
       %% Model-Based Method
        [eta, LS_MB_fail, count] = bt_LS_Relative_Shuo_MB(@(K) Func_Evaluation(K,T,A,B,Q,R,f,x0), K, grad + delta_G, Alpha, Beta, eta_min, State);


        K_plus = K - eta * (grad + delta_G);
      
        if E_x > ( 1 - Gamma ) / Gamma * norm( grad + delta_G , 'fro' )
            large_rel_error = 1;
%             disp('Large Relative Error')
        else
            large_rel_error = 0;
        end
        
        E_x_plus = ( L_r ) * norm(K_plus - K,'fro') + E_x;

        if LS_MB_fail == 0 && large_rel_error == 0
            disp('Model Base')
            new_state = 0;
            E_x = E_x_plus;
            K_list(:,:,i+1) = K_plus;
        else
            new_state = 1;
            count = 1; % For Debug use
            K_list(:,:,i+1) = K;
        end
        
        last_count = total_count;
        total_count = total_count + count;         
        MB_Start_Index(MB_count) = last_count ;
        MB_End_Index(MB_count) = total_count;
        MB_count = MB_count + 1;
        K_list_count(:,:, last_count + 1: total_count - 1 ) = reshape( repmat(K,1,total_count - last_count - 1),[p,n,total_count - last_count - 1]);
        K_list_count(:,:, total_count) = K_list(:,:,i+1);
        
       

    else
       %% Model-Free Method
       
       
       disp('Model Free')

       
       nabla_C_real = GradientEstimator(K,r,T,A,B,Q,R,f,x0);
%        keyboard
       
       last_count = total_count;
       total_count = total_count +  p * n + 1;
%         total_count = total_count  + 1
       [eta,count] = bt_line_search_MF(@(K) Func_Evaluation(K,T,A,B,Q,R,f,x0), K, nabla_C_real, Alpha, Beta);

       total_count = total_count + count;
       MF_Start_Index(MF_count) = last_count;
       MF_End_Index(MF_count) = total_count;
       
       
       K_list(:,:,i+1) = K - eta * nabla_C_real;
%        keyboard
       E_x = ( L_r ) * norm(K_list(:,:,i+1) - K,'fro');
       
       K_list_count(:,:, last_count + 1: total_count - 1 ) = reshape( repmat(K,1,total_count - last_count - 1),[p,n,total_count - last_count - 1]);
       K_list_count(:,:, total_count) = K_list(:,:,i+1);
       delta_G_list(:,:,MF_count + 1) = nabla_C_real - grad;

       MF_count = MF_count + 1;
       new_state = 0;
       
    end
    
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
    
        
end



K_list_count = K_list_count(:,:,1:total_count);

% nabla_C_list_count = [];
% 
% for i = 1:total_count
%     K = K_list_count(:,:,i);
%     nabla_C_list_count(:,:,i) = GradientEstimator(K,r,T,A,B,Q,R,f,x0);
% end
% 
% grad_norm = [];
% 
% for i = 1:total_count
%     grad_norm(i) = norm(nabla_C_list_count(:,:,i), 'fro');
% end
% 
% grad_norm = grad_norm( grad_norm > 0.001 );
% 
% keyboard

load('K_star.mat');
% keyboard

K_diff_count = K_list_count - K_star;
% K_diff = reshape(K_diff_count(:,:,1:total_count),[n*p,total_count]);


figure
for i = 1: MF_count - 1
    MF_steps = MF_Start_Index(i) - 1: MF_End_Index(i);
%     keyboard
    MF_regime = K_diff_count(:,:, MF_steps);
    MF_diff = reshape(MF_regime,[n*p,size(MF_regime,3)]);
%     keyboard
    Model_free_progeress = norms(MF_diff(:,:),[], 1);
    
    semilogy( MF_steps, Model_free_progeress , 'b','LineWidth',1.5)
    hold on
end

for i = 1: MB_count - 1
    MB_steps = MB_Start_Index(i): MB_End_Index(i);
    MB_regime = K_diff_count(:,:, MB_steps);
    MB_diff = reshape(MB_regime,[n*p,size(MB_regime,3)]);
%     keyboard
    Model_based_progeress = norms(MB_diff(:,:),[], 1);
    semilogy( MB_steps, Model_based_progeress , 'r','LineWidth',1.5)
    hold on
end


% title('Gradient Compensation Algorithm with $\eta_{min} = 5\cdot 10^{-9}$','Interpreter','latex','FontSize',16)
xlabel('Number of Function Evaluations','FontSize',16)
ylabel('$\|K - K^{*}\|$','Interpreter','latex','FontSize',16)

if eta_min == 0.05
    print('GC_LQR','-depsc2')
end

if eta_min == 0.5
    print('GC_LQR_large_eta_min','-depsc2')
end

if eta_min == 5 * 1e-6
    print('GC_LQR_small_eta_min','-depsc2')
end


if eta_min == 5 * 1e-9
    print('GC_LQR_too_small_eta_min','-depsc2')
end

end






