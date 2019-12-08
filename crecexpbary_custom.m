function [x, xs] = crecexpbary_custom(oracle, m0, x0, nu, lambda, ...
                                      curiosity_fun, tspan)
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
    baryfunc = @(t, x) baryopt(t, x, oracle, nu, lambda, curiosity_fun);
    
    x0 = [m0; x0];
    
    xhat = ode45(baryfunc, tspan, x0);
    
    % Accumulated values
    xs = xhat.y(2:end, :);
    xs = xs';
    
    % End value
    x = xs(end, :);
end

function dx = baryopt(t, x, oracle, nu, lambda, curiosity_fun)
    persistent ei_s z_s xhats ms;
    
    if(isempty(ei_s))
        ei_s = [];
        z_s = [];
        xhats = [];
        ms = [];
    end
    
    m  = x(1);
    xhat = x(2:end);
    
    x = xhat + curiosity_fun(t, xhat, m);
    e_i = exp(-nu*oracle(x));
    
    dm = exp(-lambda^t)*e_i;
    
    dxhat = (1/m)*(x - xhat*exp(-lambda^t))*e_i;

    dx = [dm; dxhat];
    
    ei_s = [ei_s; e_i];
    xhats = [xhats, xhat];
    ms = [ms; m];
    
    assignin('base', 'exps', ei_s)
    assignin('base', 'xhats', xhats)
    assignin('base', 'ms', ms)
end