function [x, xs, m] = d2recexpbary_custom(oracle, m0, x0, nu, sigma, ...
                                         lambda, iterations)
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
    
    oracle = @(x) dsecond_order_fun(x, oracle);

    % Initialization    
    xhat_2 = x0;
    m_2 = m0;    
    z_bar_2 = zeros(size(x0));
    
    delta_x_2 = zeros(size(x0));
    z_bar_1 = zeros(size(x0));
    
    z_2 = normrnd(z_bar_2, sigma);
    x_2 = xhat_2 + z_2;
    f_2 = oracle(x_2);
    
    % Previous value of function
    z_1 = normrnd(z_bar_1, sigma);
    x_1 = xhat_2 + z_1;
    
    f_1 = oracle(x_1);
    e_1 = exp(-nu*f_1);
    
    % Mass components    
    m_1 = m_2 + e_1;    
    xhat_1 = (1/m_1)*(lambda*m_2*xhat_2 + x_1*e_1);
    delta_x_1 = xhat_1 - xhat_2;
    
    xs = [];
    ms = [];
    
    solution_found = false;
    
    i = 1;
    while(~solution_found)
       % Calculation of mean for stochastic function
       z_bar = (1/m_1)*(exp(-nu*f_1)*delta_x_1 + z_bar_1*m_2);
       z = normrnd(z_bar, sigma);
        
       % Current value of positoin
       x = xhat_1 + z;
       
       % Mass component
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
        
       % Barycenter point
       xhat = (1/m)*(lambda*m_1*xhat_1 + x*e_i);
        
       solution_found = i >= iterations;
       
       % Previous elements
       deltax_1 = xhat - xhat_1;
       
       % Updates
       m_2 = m_1;
       m_1 = m;
       
       xhat_2 = xhat_1;
       xhat_1 = xhat;
       
       xs = [xs; xhat'];
       ms = [ms; m];
    
        i = i + 1;
    end
    
    xs = [x0'; xs];
end
