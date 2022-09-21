function [t,State,count] = bt_LS_MB(f, x, grad, rho, Beta, eta_min, State)
t = 1;
dx = -grad;



% (TODO): still, x+t*dx may go outside the domain of f

% Backtracking line search

count = 2;

while f(x + t*dx) >= f(x) + rho * t * grad(:)' * dx(:)
    
    if  t >= eta_min
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


