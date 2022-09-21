function x_star = ComputeOptimal(f, hat_x_star)

n = length(hat_x_star);
cvx_begin
    cvx_precision high
    variable x(n)
    minimize( f(x) )
cvx_end
x_star = x; 
% keyboard
% save x_star.mat
end

