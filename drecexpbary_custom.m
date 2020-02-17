function [x, xs, m] = drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
                                         lambda, lambda_z, iterations, ...
                                         is_accel)
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
    m_2 = 0;
    m_1 = m0;
    xhat_1 = x0;
    delta_x_1 = zeros(size(x0));
    z_bar_1 = zeros(size(x0));
    f_1 = oracle(x0);
    
    xs = [];
    ms = [];
    
    solution_found = false;
    
    i = 1;
    zbars = [];
    while(~solution_found)
       % Calculation of mean for stochastic function
       e_z = exp(-nu*f_1);
       
       if(is_accel)
           z_bar = (1/m_1)*(-lambda_z*e_z*delta_x_1 + z_bar_1*m_2);
       else
           z_bar = zeros(size(delta_x_1));
       end
       
       z_bar = vpa(z_bar);
       zbars = [zbars; z_bar'];
       assignin('base', 'zbars', zbars);
       
       z = normrnd(double(z_bar), double(sigma));
        
       % Current value of positoin
       x = xhat_1 + z;
       
       % Mass component
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
        
       % Barycenter point
       xhat = (1/m)*(lambda*m_1*xhat_1 + x*e_i);
        
       solution_found = i >= iterations;
       
       % Previous elements
       delta_x_1 = xhat - xhat_1;
       
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
