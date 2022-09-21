function [hat_nabla_C,hat_C] = GradientEstimator(K,r0,T,A,B,Q,R,f,x0)

% if abs(eigs(A - B*K, 1)) >= 1
%     hat_C = Inf;
%     keyboard
% else

    [p,n] = size(K);
    
    np = n*p; M = 2*np;
    r = r0 * reshape([eye(np), -eye(np)], [p n M]);
    rad = norm( r(:,:,1) );
    CjUj = zeros(p,n);
    Cj = 0;
    for k = 1:n
        for i = 1:M
        x = x0(:,k);
        %% Compute hat_C( K + r )
        hat_C = 0;
            for j = 1:T
                u = - ( K + r(:,:,i) ) * x;
                hat_C = hat_C + x'*Q*x + u'*R*u;
                x_next = A*x + B*u + f(x);
                x = x_next;
            end
        %%
        Cj_next  = Cj + hat_C;
        CjUj_next = CjUj + hat_C * r(:,:,i);

        Cj = Cj_next;
        CjUj = CjUj_next;

        end
    end

    
    hat_C = 1/(n*M) * Cj;
    hat_nabla_C = 1/(M*n) * np / (rad^2) * CjUj;
end