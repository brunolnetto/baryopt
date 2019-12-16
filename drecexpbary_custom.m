function [x, xs] = drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
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
    
    m_1 = m0;
    xhat_1 = x0;
    
    xs = [];
    ms = [];
    
    solution_found = false;
    
    i = 1;
    while(~solution_found)
       z0 = zeros(size(xhat_1));
       z = normrnd(z0, sigma);
       
       x = xhat_1 + z;
       e_i = exp(-nu*oracle(x));
       m = lambda*m_1 + e_i;
       xhat = (1/m)*(lambda*m_1*xhat_1 + x*e_i);
        
       solution_found = i >= iterations;
       
       % Updates
       m_1 = m;
       deltax_1 = xhat - xhat_1;
       xhat_1 = xhat;
       
       xs = [xs; xhat'];
       ms = [ms; m];
    
        i = i + 1;
    end
    
    xs = [x0'; xs];
end
