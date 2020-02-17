function [x, xs, zs] = crecexpbary_custom(oracle, m0, x0, zbar_0, ...
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
    baryfunc = @(t, x) baryopt(t, x, oracle, nu, lambda, lambda_z, sigma);
    
    if(m0 == 0)
        error('m0 MUST be different from 0!');
    end
    
    n = length(x0);
    x0 = [m0; m0; x0; zbar_0];
    
    xhat = my_ode45(baryfunc, tspan, x0);
    
    % Accumulated values
    xs = xhat(3:2+n, :);
    zs = xhat(3+n:end, :);
    
    xs = xs';
    zs = zs';
    
    % End value
    x = xs(end, :);
end

function dx = baryopt(t, x, oracle, nu, lambda, lambda_z, sigma)
    persistent xhats e_s ms counter;
    
    if(isempty(counter))
        counter = 0;
        xhats = [];
        ms = [];        
    end
        
    n = (length(x)-2)/2;
    m  = x(1);
    m_z  = x(2);
    xhat = x(3:2+n);
    % zbar = x(3+n:end);
    zbar = zeros(size(xhat));
        
    x = normrnd(zbar, sigma);
    e_i = exp(-nu*oracle(x));
    
    dm = exp(-lambda^t)*e_i;
    dm_z = exp(-lambda_z^t)*e_i;
    dxhat = (1/m)*(x - xhat)*e_i;
    % dzbar = (1/m_z)*(dxhat - zbar)*exp(-lambda_z^t)*e_i;
    dzbar = zeros(size(xhat));
    
    dx = [dm; dm_z; dxhat; dzbar];
    
    counter = counter + 1;
    
    if(counter == 1)
        xhats = [xhats, xhat];
        ms = [ms; m];
        e_s = [e_s; e_i];
                
        assignin('base', 'xhats', xhats);
        assignin('base', 'ms', ms);
        assignin('base', 'e_s', e_s);
    elseif(counter == 4)
        counter = 0;
    end
end