function P = rand_psdmat(n)

% Return a random (symmetric) positive definite matrix

P1 = randn(n);
P = P1'*P1;

end