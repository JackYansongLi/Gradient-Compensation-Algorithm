function grad = Model_Gradient_Estimator(K,A,B,Q,R,sigma,mu2)
    Pk = dlyap( (A-B*K)' , Q + K'*R*K );
    Ek = (R+B'*Pk*B)*K-B'*Pk*A;
    sigma_k = dlyap((A-B*K),sigma+mu2);
    grad = 2 * Ek * sigma_k;
end