function [x_star, xhats, xs, ms, deltas, zbars, Fs_n, Fbars_n] = ...
            drecexpbary_custom(oracle, m0, xhat0, nu, sigma, ...
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
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    xhat_1 = zeros(size(xhat0));
    
    xhats = xhat0';
    ms = m0;
    deltas = [];
    Fs_n = [];
    Fbars_n = [];
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
       x = xhat_1 + z;
       
       
       % Mass component
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
       
       F_n = vpa(e_i/(m_1 + e_i));
       Fbar_n = (m_1/(m_1 + e_i))*F_n;
       
       % Barycenter point
       xhat = vpa(xhat_1 + F_n*z);              
       solution_found = i >= iterations;
       
       % Updates
       m_1 = m;
       delta_xhat_1 = xhat - xhat_1;
       xhat_1 = xhat;
       
       xs = [xs; x'];
       xhats = [xhats; xhat'];
       ms = [ms; m];
       deltas = [deltas; delta_xhat_1'];
       Fs_n = [Fs_n; F_n];
       Fbars_n = [Fbars_n; Fbar_n];
       
       i = i + 1;
    end
    
    x_star = xhats(end, :);
end
