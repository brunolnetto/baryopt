function z = curiosity_fun(t, xhat, sigma)
    z = normrnd(zeros(size(xhat)), sigma);
end