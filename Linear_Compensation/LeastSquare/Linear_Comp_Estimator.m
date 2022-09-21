function [W,b] = Linear_Comp_Estimator(K,r0,model_gradient,true_gradient)

    [p,n] = size(K);
    np = n*p; M = 2*np;
%     keyboard
    r = r0 * reshape([eye(np), -eye(np)], [p n M]);
    
    %%Compute Matrix X
    K0 = [1;reshape(K, [np 1])];
    K_reshape = reshape(K + r, [p*n M]);
    X_extend = [ones(M,1), K_reshape'];
    X_extend = [K0'; X_extend];
    
    %%Compute Matrix Y
    Y = zeros(M,np);
    for i = 1:M
       Ki = K + r(:,:,i);
       yi = true_gradient(Ki) - model_gradient(Ki);
       Y(i,:) = reshape(yi, [1, np]);
    end
    y0 = true_gradient(K) - model_gradient(K);
    y0_reshape = reshape(y0, [1, np]);
    Y = [y0_reshape ; Y];
    
    W_extend = Y'*X_extend/(X_extend'*X_extend);
    
    W = W_extend(:,2:end);
    b = W_extend(:,1);
    
%     W = 0;
%     b = GradientEstimator(K,r0,T,A,B,Q,R,f,x0) - Model_Gradient_Estimator(K,A,B,Q,R,sigma,mu2);
%     b = reshape(b, [np,1]);

end