function [x, xs, zs] = crecexpbary_custom(oracle, m0, x0, zbar_0, ...
                                          nu, lambda, sigma, tspan)
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
    baryfunc = @(t, x) baryopt(t, x, oracle, nu, lambda, sigma);
    
    if(m0 == 0)
        error('m0 MUST be different from 0!');
    end
    
    n = length(x0);
    x0 = [m0; x0; zbar_0];
    
    xhat = my_ode45(baryfunc, tspan, x0);
    
    % Accumulated values
    xs = xhat(2:1+n, :);
    zs = xhat(2+n:end, :);
    
    xs = xs';
    zs = zs';
    
    % End value
    x = xs(end, :);
end

function dx = baryopt(t, x, oracle, nu, lambda, sigma)
    persistent ei_s z_s xhats ms;
    
    if(isempty(ei_s))
        ei_s = [];
        z_s = [];
        xhats = [];
        ms = [];
    end
    
    n = length((length(x)-1)/2);
    
    m  = x(1);
    xhat = x(2:2+n);
    zbar = x(3+n:end);
    
    x = normrnd(zbar, sigma);
    e_i = exp(-nu*oracle(x));
    
    dm = exp(-lambda^t)*e_i;
    dxhat = (1/m)*(x - xhat*exp(-lambda^t))*e_i;
    dzbar = (dxhat - zbar)*e_i/m;
%     dzbar = zeros(size(xhat));
    
    dx = [dm; dxhat; dzbar];
    
    ei_s = [ei_s; e_i];
    xhats = [xhats, xhat];
    ms = [ms; m];
    
    assignin('base', 'exps', ei_s)
    assignin('base', 'xhats', xhats)
    assignin('base', 'ms', ms)
end