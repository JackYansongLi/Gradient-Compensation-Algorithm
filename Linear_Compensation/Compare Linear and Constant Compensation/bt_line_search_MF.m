function [t,count] = bt_line_search_MF(f, x, grad, rho, Beta)
t = 2;
dx = -grad;

% Choose t such that x + t*dx stays within the domain of f
while x + t*dx < 0
    t = Beta * t;
end

% (TODO): still, x+t*dx may go outside the domain of f

% Backtracking line search

count = 2;

while f(x + t*dx) - f(x) >=   rho * t * grad(:)' * dx(:)
    if t > 1e-10
        t = Beta * t;
        count = count + 1;
    else
        disp('Line Search Fail for Model-Free')
        keyboard;
        break;
    end
end


end