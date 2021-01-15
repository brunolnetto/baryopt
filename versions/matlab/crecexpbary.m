function [x, xs] = crecexpbary(oracle, m0, x0, nu, Sigma, ...
                               lambda, zeta, tspan)
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
    baryfunc = @(t, x) baryopt(t, x, oracle, nu, Sigma, lambda, zeta);
    
    x0 = [m0; x0];
    
    xhat = ode45(baryfunc, tspan, x0);
    
    % Accumulated values
    xs = xhat.y(2:end, :);
    xs = xs';
    
    % End value
    x = xs(end, :);
end

function dx = baryopt(t, x, oracle, nu, Sigma, lambda, zeta)
    persistent ei_s z_s xhats ms;
    
    if(isempty(ei_s))
        ei_s = [];
        z_s = [];
        xhats = [];
        ms = [];
    end
    
    m  = x(1);
    xhat = x(2:end);
    
    zhat = exp(-nu*xhat);
    z = normrnd(zhat, Sigma);
    
    x = xhat + z;
    e_i = exp(-nu*oracle(x));
    
    dm = (lambda^t)*e_i;
    
    dxhat = (1/m)*(x - xhat*exp(-lambda^t))*e_i;

    dx = [dm; dxhat];
    
    ei_s = [ei_s; e_i];
    z_s = [z_s; z];
    xhats = [xhats, xhat];
    ms = [ms; m];
    
    assignin('base', 'exps', ei_s)
    assignin('base', 'z_s', z_s)
    assignin('base', 'xhats', xhats)
    assignin('base', 'ms', ms)
end