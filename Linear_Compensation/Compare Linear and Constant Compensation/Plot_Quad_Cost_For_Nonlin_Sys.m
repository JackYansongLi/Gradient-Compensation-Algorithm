function Plot_Quad_Cost_For_Nonlin_Sys
    
    rng(60607);
    T = 50;                  % Factor of Function Evaluation, i.e., hat_C
    n = 1; p = 1; np = n*p;
    Q = 2*eye(n); R = eye(p);
    A = randn(n,n);  B = randn(n,p); N = zeros(n,p);
    [K_lin] = dlqr(A,B,Q,R,N);
    mf_factor = 0.01;
    f = @(x) (mf_factor*x/(1-0.9*sin(x)));
    x0 = randn(n,n);
    
    x_list = K_lin + linspace(-1,1,100);
    y_list = zeros(size(x_list));
    for i = 1:length(x_list)
        y_list(i) = Func_Evaluation(x_list(i),T,A,B,Q,R,f,x0);
    end
    
    plot(x_list,y_list)



end