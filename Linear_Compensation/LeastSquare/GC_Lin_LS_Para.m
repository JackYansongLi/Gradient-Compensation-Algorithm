rng(60607)
n = 4;
c1 = 0.01; %Nonlinear Factor
c2 = 0.01; %Drifting Factor
A = randn(n,n); 
b = randn(n,1);
f = @(x)( pow_pos(norm(A*x-b),2) + c1 * pow_pos(norm( A * (x-c2) - b ),2) );
d_f = @(x) ( 2 * A' * ( A * x - b ) + 2 * c1 * A' * ( A * ( x - c2 ) - b ) );
hat_f = @(x) ( pow_pos(norm( A * x - b ),2) );
d_hat_f = @(x) ( 2 * A' * ( A * x - b ) );
hat_x_star = (A'*A)\(A'*b);

x0 = ones(n,1) ; % Ini_value be the optimal solution of the model part


x_star = ComputeOptimal(f, hat_x_star);



x_list = [];
x_list(:,1) = x0;
total_count = 1;
x_list_count = zeros(n,1000);
x_list_count(:,1) = x0;
r0 = 0.001; % Factor of Linear Compensation Estimator
rho = 0.3; Beta = 0.5;  % Factor of Line search
MF_count = 1;
MB_count = 1;
new_state = 0;
MF_Start_Index = [];
MF_End_Index = [];
MB_Start_Index = [];
MB_End_Index = [];
W_list = [];
b_list = [];
W_list(:,:,1) = zeros(n,n);
b_list(:,1) = zeros(n,1);
tol = 1e-8;
i = 0;
x = x0; 