function A = rand_condmat(n,c)

% Return a random (symmetric) positive definite matrix with condition
% number c

a = rand_psdmat(n);

[U,S,V] = svd(a);
S(S~=0) = linspace(c,1,n);

A = U*S*V';

end