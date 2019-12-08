function z = curiosity_fun(t, xhat, m, oracle, nu, sigma, zeta)
    z0 = zeros(size(xhat));

    % Accelerates the convergence
    acceleration_curiosity = normrnd(z0, sigma);
    x = xhat + acceleration_curiosity;

    e_i = exp(-nu*oracle(x));
    
    xhat_p = (1/m)*(x - xhat)*e_i;
    
    zbar = z0;
%     zbar = -zeta*xhat_p;
    
    z = normrnd(zbar, sigma);
end