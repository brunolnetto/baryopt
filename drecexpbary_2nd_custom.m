function [x_star, xhats, xs, ...
          ys, ms, deltas, zbars] = ...
            drecexpbary_2nd_custom(oracle, m0, xhat0, nu, sigma, ...
                                   lambda, iterations, accel_fun)
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
    
    PRECISION = 10;
    digits(PRECISION);
    
    % Mass components
    m_1 = m0;
    xhat_2 = xhat0;
    xhat_1 = xhat0;
    y_1 = xhat0;
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    
    xhats = xhat0';
    ys = xhat0';
    ms = m0;
    deltas = [];
    ys = [];
    xs = [];
    
    solution_found = false;
    
    i = 1;
    zbars = [];
    while(~solution_found)
       
       % Calculation of mean for stochastic function
       zbar = accel_fun(m_1, xhat_1, delta_xhat_1);
       
       zbars = [zbars; zbar'];
       z = normrnd(double(zbar), double(sigma));
       
       % Current value of position
       x = y_1 + z;       
       
       % Mass component
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
       
       % Barycenter point
       xhat = (1/m)*(m_1*xhat_1 + e_i*x);
       
       theta_i = (i - 1)/(i + 2);
       y = xhat_1 + theta_i*(xhat_1 - xhat_2);
       
       solution_found = i >= iterations;
       
       % Updates
       m_1 = m;
       delta_xhat_1 = xhat - xhat_1;
       xhat_2 = xhat_1;
       xhat_1 = xhat;
       y_1 = y;
       
       xs = [xs; x'];
       xhats = [xhats; xhat'];
       ys = [ys; y'];
       ms = [ms; m];
       deltas = [deltas; delta_xhat_1'];
       
       i = i + 1;
    end
    
    x_star = xhats(end, :);
end
