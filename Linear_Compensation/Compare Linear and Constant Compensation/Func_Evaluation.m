function hat_C = Func_Evaluation(K,T,A,B,Q,R,f,x0)

    hat_C = 0;
    n = size(x0(:,1),1);
    
    for i = 1:n
        x = x0(:,i);
        for j = 1:T
            u = - K  * x;
            hat_C = hat_C + x'*Q*x + u'*R*u;

            x_next = A*x + B*u + f(x);

            x = x_next;
        end
    end

    hat_C = 1/n * hat_C;

end
