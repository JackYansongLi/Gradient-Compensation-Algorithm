function Qu_LQR

rng(60607)
iters = 850; 

Alpha = 0.3; Beta = 0.5; % Factor of Line Search

T = 50;                  % Factor of Function Evaluation, i.e., hat_C

r = 0.001;               % Gradient Estimation radius


n = 4;
p = 3;
Q = 2*eye(n);
R = eye(p);


A = randn(n,n); 
B = randn(n,p);
N = zeros(n,p);

[K_lin] = dlqr(A,B,Q,R,N);

mf_factor = 0.01;

f = @(x) [mf_factor*x(1)/(1-0.9*sin(x(1))); mf_factor*x(2)/(1-0.9*sin(x(2))); ...
          mf_factor*x(3)/(1-0.9*sin(x(3))); mf_factor*x(4)/(1-0.9*sin(x(4)))];



total_count = 0;
K_list = zeros(p,n,iters);
K_list(:,:,1) = K_lin;
K_list_count = zeros(p,n,1000);
K_list_count(:,:,1) = K_lin;

x0 = randn(n,n);
K = K_lin;
load('K_star.mat')
i = 0;
tol = 1e-4;
% for i = 1:iters-1
while norm(K - K_star, 'fro') > tol
    i = i+1;
    
    K = K_list(:,:,i);
    
    
    grad = GradientEstimator(K,r,T,A,B,Q,R,f,x0);
%     Func_Evaluation(K,T,A,B,Q,R,f,x0);    
%     keyboard
    
    last_count = total_count;
    total_count = total_count + p * n + 1;
    [eta,count] = bt_line_search_MF(@(K) Func_Evaluation(K,T,A,B,Q,R,f,x0), K, grad, Alpha, Beta);

    total_count = total_count + count;
    K_list(:,:,i+1) = K - eta*grad;
    K_list_count(:,:, last_count + 1: total_count - 1 ) = reshape( repmat(K,1,total_count - last_count - 1),[p,n,total_count - last_count - 1]);
    K_list_count(:,:, total_count) = K_list(:,:,i+1);
    
    if mod(i,5) == 0
        disp(['Iteration: ',num2str(i)])
    end
end



K_list_count = K_list_count(:,:,1:total_count);
% K_star = K_list_count(:,:,total_count);
% keyboard

K_diff_count = K_list_count - K_star;
K_diff_reshape = reshape(K_diff_count(:,:,1:total_count),[n*p,total_count]);

vert_axis = norms(K_diff_reshape ,[], 1);

figure
semilogy( 1:size(vert_axis,2), vert_axis,'LineWidth',1.5)
% title('Qu et al.: Modified LQR','FontSize',16)
xlabel('Number of Function Evaluations','FontSize',16)
ylabel('$\|K - K^{*}\|$','Interpreter','latex','FontSize',16)

print('Qu_LQR','-depsc2')
end
% keyboard;






