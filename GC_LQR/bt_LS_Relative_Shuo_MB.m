function [t,State,count] = bt_LS_Relative_Shuo_MB(f, x, grad, Alpha, Beta, eta_min, State)
dx = -grad;

% Choose t such that x + t*dx stays within the domain of f
% while abs(eigs(A - B*(x + t*dx), 1)) >= 1
%     t = Beta * t;
% end

% (TODO): still, x+t*dx may go outside the domain of f

% Backtracking line search

count = 1;
t = eta_min;

if f(x + t*dx) >= f(x) + Alpha * t * grad(:)' * dx(:)
    State = 1;
%     disp('Line Search Failed')
    return;
else
    t = 1;
end


while f(x + t*dx) >= f(x) + Alpha * t * grad(:)' * dx(:)
    
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


