function [x, xs, m, deltas, zbars] = ...
            drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
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
        
    % Mass components
    m_1 = m0;
    z_bar_1 = zeros(size(x0));
    delta_xhat_1 = zeros(size(x0));
    xhat_1 = zeros(size(x0));
    
    xs = [];
    ms = [];
    deltas = [];
    
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
       
       % Barycenter point
       sum_hat_1 = m_1*xhat_1;
       xhat = (1/m)*(lambda*sum_hat_1 + x*e_i);
        
       solution_found = i >= iterations;
       
       % Updates
       m_1 = m;
       delta_xhat_1 = xhat - xhat_1;
       xhat_1 = xhat;
       
       xs = [xs; xhat'];
       ms = [ms; m];
       deltas = [deltas; delta_xhat_1']; 
       
       i = i + 1;
    end
    
    xs = [x0'; xs];
end
