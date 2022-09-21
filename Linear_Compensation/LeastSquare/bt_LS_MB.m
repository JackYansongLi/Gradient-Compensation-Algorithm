function [t,State,count] = bt_LS_MB(f, x, grad, rho, Beta, State)
t = 1;
dx = -grad;

% Choose t such that x + t*dx stays within the domain of f
while x + t*dx <= 0
    t = Beta * t;
end

% (TODO): still, x+t*dx may go outside the domain of f

% Backtracking line search

count = 2;

while f(x + t*dx) >= f(x) + rho * t * grad(:)' * dx(:)
    
    if  t >= 0.01
        t = Beta * t;
        count = count + 1;
        State = 0;
    else
        State = 1;
        disp('Line Search Failed')
%         keyboard;
        break;
    end

end

end


