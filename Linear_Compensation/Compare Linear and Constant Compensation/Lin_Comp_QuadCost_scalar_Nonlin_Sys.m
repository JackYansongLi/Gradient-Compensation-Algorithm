function Lin_Comp_QuadCost_scalar_Nonlin_Sys
iters = 300;
rng(60607);
Alpha = 0.3; Beta = 0.5; eta_min = 5*1e-4;  % Factor of Line Search
r0 = 0.001;               % Gradient Estimation radius
T = 1;                  % Factor of Function Evaluation, i.e., hat_C
n = 1; p = 1; np = n*p;
Q = 2*eye(n); R = eye(p);
A = randn(n,n);  B = randn(n,p); N = zeros(n,p);
K_lin = dlqr(A,B,Q,R,N);
mf_factor = 0.01;
f = @(x) (mf_factor*x/(1-0.9*sin(x)));
x0 = randn(n,n);
mu = mean(x0,2);
mu2 = mu*mu';
sigma = cov(x0,1); 

MF_count = 1;
MB_count = 1;
new_state = 0;
total_count = 1;
W_list = [];
b_list = [];
W_list(:,:,1) = zeros(np,np);
b_list(:,1) = zeros(np,1);

K_list = zeros(p,n,iters);
K_list(:,:,1) = K_lin;
K_list_count = zeros(p,n,1000);
K_list_count(:,:,1) = K_lin;

MF_Start_Index = [];
MF_End_Index = [];
MB_Start_Index = [];
MB_End_Index = [];

for i = 1:iters-1
    K = K_list(:,:,i);
    
    W = W_list(:,:,MF_count);
    b = b_list(:,:,MF_count);
    K_reshape = reshape(K, [np, 1]);
    delta_G_reshape = W*K_reshape + b;
    delta_G = reshape(delta_G_reshape, [p, n]);
%     keyboard
    
    grad = Model_Gradient_Estimator(K,A,B,Q,R,sigma,mu2);

    
    
    State = new_state;

    if State == 0
       %% Model-Based Method
        
        [eta, LS_MB_fail, count] = bt_LS_MB(@(K) Func_Evaluation(K,T,A,B,Q,R,f,x0), K, grad + delta_G, Alpha, Beta, eta_min, State);

        K_plus = K - eta * (grad + delta_G);

        if LS_MB_fail == 0 
            disp('Model Base')
            new_state = 0;
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
       nabla_C_real = GradientEstimator(K,r0,T,A,B,Q,R,f,x0);
       
       last_count = total_count;
       total_count = total_count +  p * n + 1;
       [eta,count] = bt_line_search_MF(@(K) Func_Evaluation(K,T,A,B,Q,R,f,x0), K, nabla_C_real, Alpha, Beta);
       total_count = total_count + count;
       MF_Start_Index(MF_count) = last_count;
       MF_End_Index(MF_count) = total_count;
       
       
       K_list(:,:,i+1) = K - eta * nabla_C_real;
       
       K_list_count(:,:, last_count + 1: total_count - 1 ) = reshape( repmat(K,1,total_count - last_count - 1),[p,n,total_count - last_count - 1]);
       K_list_count(:,:, total_count) = K_list(:,:,i+1);
       [W,b] = Linear_Compensation_Estimator(K,r0,T,A,B,Q,R,f,x0,sigma,mu2);
       W_list(:,:,MF_count + 1) = W;
       b_list(:,:,MF_count + 1) = b;
       MF_count = MF_count + 1;
       new_state = 0;
       
    end
    
    
    if mod(i,10) == 0
        disp(['Iteration: ',num2str(i)])
    end    
end




K_list_count = K_list_count(:,:,1:total_count);
K_star =  K_list_count(:,:,total_count);
K_diff_count = K_list_count - K_star;

figure
for i = 1: MF_count - 1
    MF_steps = MF_Start_Index(i) - 1: MF_End_Index(i);
    MF_regime = K_diff_count(:,:, MF_steps);
    MF_diff = reshape(MF_regime,[n*p,size(MF_regime,3)]);
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



end