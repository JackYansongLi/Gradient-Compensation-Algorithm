function [t,count] = bt_line_search_MF(f, x, grad, rho, Beta)
% t = 10;
t = 1;
dx = -grad;


% Choose t such that x + t*dx stays within the domain of f
% while abs(eigs(A - B*(x + t*dx), 1)) >= 1
%     t = Beta * t;
% end

% (TODO): still, x+t*dx may go outside the domain of f

% Backtracking line search

count = 1;


% Define a function handle delta_f

while f(x + t*dx) - f(x) -  rho * t * grad(:)' * dx(:) >= 0
        t = Beta * t;
        count = count + 1;
end

end