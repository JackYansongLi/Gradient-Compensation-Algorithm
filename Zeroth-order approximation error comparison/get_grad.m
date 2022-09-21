function [g]=get_grad(f, x, delta)
    % f The function we want to compute the gradient
    % x the point where we want to compute
    % n the size of the square matrix
    % M sample size
if any(isnan(x))
    warning('x is NaN!');
    keyboard; 
end
    [n,p] = size(x);
    np = n*p; M = 2*np;
    Delta = delta * reshape([eye(np), -eye(np)], [n p M]);
    xi = x + Delta;
    rx = reshape(x, [np,1]);

    cvx_begin quiet
        variable grad(np) 
        expression r(M)
        for i = 1:M
            rxi = reshape(xi(:,:,i), [np,1]);
            r(i) = f(xi(:,:,i)) - (f(x) + grad'*(rxi - rx));
        end
        minimize (norm(r))
    cvx_end
    g = grad;
end