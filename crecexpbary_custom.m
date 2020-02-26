function [x, xs, zs, ms] = crecexpbary_custom(oracle, m0, x0, zbar_0, ...
                                          nu, lambda, lambda_z, ...
                                          sigma, tspan)
% Recursive barycenter algorithm for direct optimization
% https://arxiv.org/abs/1801.10533
% In:
%   - oracle [function]: Oracle function
%   - x0 []: Query values
%   - nu []: positive value (Caution on its value due overflow)
%   - Sigma []: Std deviation of normal distribution
%   - lambda []: Forgetting factor
%   - time []: Time lapse
% Out:
%   - x []: Optimum position
%   - xs []: Optimum position evolution
    degree = 8;
    baryfunc = @(t, x) baryopt(t, x, oracle, ...
                               nu, lambda, lambda_z, sigma, degree);
    
    if(m0 == 0)
        error('m0 MUST be different from 0!');
    end
    
    n = length(x0);
    x0 = [m0; m0; x0; zbar_0];
    
    [tspan, sol] = ode(degree, baryfunc, x0, tspan);
    
    % Accumulated values
    ms = sol(1, :);
    xs = sol(3:2+n, :);
    zs = sol(3+n:end, :);
    
    xs = xs';
    zs = zs';
    
    % End value
    x = xs(end, :);
end

function dx = baryopt(t, x, oracle, nu, lambda, lambda_z, sigma, degree)
    persistent counter t_hat dxhats;
    
    if(isempty(counter))
        counter = 0;
        t_hat = [];
        dxhats = [];
    end
    
    n = (length(x)-2)/2;
    m  = x(1);
    m_z  = x(2);
    
    xhat = x(3:2+n);
    zbar = zeros(size(xhat));
        
    z = normrnd(zbar, sigma);
    x = xhat + z;
    
    e_i = exp(-nu*oracle(x));
    
    dm = exp(-lambda^t)*e_i;
    dm_z = exp(-lambda_z^t)*e_i;
    
    dxhat = (1/m)*(x - xhat)*e_i;
    dzbar = (1/m_z)*(dxhat - zbar)*exp(-lambda_z^t)*e_i;
    
    dx = [dm; dm_z; dxhat; dzbar];

    counter = counter + 1;
    if(counter == 1)
        dxhats = [dxhats, dx];
        t_hat = [t_hat, t];
        assignin('base', 'dxhats', dxhats);
        assignin('base', 't_hat', t_hat);
    end

    if(counter == degree)
        counter = 0;
    end
end