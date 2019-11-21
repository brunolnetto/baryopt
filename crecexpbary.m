function [x, xs] = crecexpbary(oracle, x0, nu, Sigma, lambda, time)
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
    baryfunc = @(t, x) cbaryopt(t, x, oracle, nu, Sigma, lambda);
    
    x0 = [0; x0];
    
    xhat = ode45(baryfunc, time, x0);
    xs = xhat(:, end);
    x = xhat(end, :);
end

function dx = cbaryopt(t, x, oracle, nu, Sigma, lambda)
    z = normrnd(0, Sigma);
    
    persistent ei_s z_s xhats ms;
    
    if(isempty(ei_s))
        ei_s = [];
        z_s = [];
        xhats = [];
        ms = [];
    end
    
    m  = x(1);
    xhat = x(2:end);
    
    x = xhat + z;
    e_i = exp(-nu*oracle(x));
    dm = lambda*e_i;
    dxhat = (1/m)*(x - lambda*xhat)*e_i;
    
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