function [x_star, xhats, xs, ms, ...
          deltas, zbars, Fs_n, Fbars_n] = ...
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
    m_1 = terop(m0 ~= 0, m0, eps);
    xhat_1 = xhat0;
    x_1 = xhat0;
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    
    xhats = xhat0';
    ms = m0;
    deltas = [];
    zbars = zeros(size(xhat0'));
    Fs_n = [];
    Fbars_n = [];
    xs = [];
    
    solution_found = false;
    
    i = 1;
    while(~solution_found)
       
       % Calculation of mean for stochastic function
       zbar = accel_fun(m_1, x_1, delta_xhat_1);
       
       z = zeros(length(zbar), 1);
       for j = 1:length(zbar)
           z(j) = gaussianrnd(double(zbar(j)), ...
                              double(sigma));
       end
       
       % Current value of position
       x = xhat_1 + z;       
       
       % Mass component
       
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
       
       F_n = vpa(e_i/m);
       Fbar_n = (m_1/m)*F_n;
       
       % Barycenter point
       xhat = vpa(xhat_1 + F_n*z);              
       
       solution_found = i >= iterations;
       
       % Updates
       m_1 = m;
       delta_xhat_1 = xhat - xhat_1;
       xhat_1 = xhat;
       x_1 = x;
       
       xs = [xs; x'];
       xhats = [xhats; xhat'];
       ms = [ms; m];
       deltas = [deltas; delta_xhat_1'];
       Fs_n = [Fs_n; F_n];
       Fbars_n = [Fbars_n; Fbar_n];
       zbars = [zbars; zbar'];
       
       i = i + 1;
       
       disp(sprintf('i=%3d: fval: %.3f', i, oracle(x)));
    end
    
    x_star = xhats(end, :);
end
