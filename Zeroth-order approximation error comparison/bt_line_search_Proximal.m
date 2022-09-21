function [t,count] = bt_line_search_Proximal(g , x, grad, Beta, P, n)


% f(x) = g(x) + h(x)
% \nabla g(x) = grad
% x = prox_h ( x - t * grad )

t = 1;
dx = grad;


u  = x - t * grad;

prox = ( P + 1/t * eye(n) ) \ ( 1/t * u );



Gt = 1/t * (x- prox);

while g(x - t*Gt) >= g(x) -  t * grad(:)' * Gt + t/2 * (Gt') * Gt

    t = Beta * t;
    
    u  = x - t * grad;

    prox = ( P + 1/t * eye(n) ) \ ( 1/t * u );

    Gt = 1/t * (x- prox);

end


end