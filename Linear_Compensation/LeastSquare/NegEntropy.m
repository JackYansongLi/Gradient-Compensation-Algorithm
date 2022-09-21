function neg_entropy = NegEntropy(x)
    neg_entropy = vpa(x'*log(x));
end