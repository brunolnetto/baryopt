function [x, xs] = drecexpbary(oracle, x0, nu, sigma, zeta, lambda, ...
                              iterations)
% Recursive barycenter algorithm for direct optimization
% In:
%   - oracle [function]: Oracle function
%   - x0 []: Initial query values
%   - nu []: positive value (Caution on its value due overflow)
%   - sigma []: Std deviation of normal distribution
%   - zeta []: Proportional value for mean of normal distribution
%   - lambda []: Forgetting factor
%   - iterations []: Maximum number of iterations
% Out:
%   - x []: Optimum position
%   - xs []: Optimum position evolution
    
    xhat_1 = x0;
    m_1 = 0;
    deltax_1 = zeros(size(x0));
    solution_found = false;    
    
    fis = oracle(xhat_1);
    xs = [];
    
    i = 1;
    while(~solution_found)
        i = i + 1;
        z = normrnd(zeta*deltax_1, sigma);
        
        x = xhat_1 + z;
        
        e_i = exp(-nu*oracle(x));
        
        m = lambda*m_1 + e_i;
        xhat = (1/m)*(lambda*m_1*xhat_1 + x*e_i);
        
        xs = [xs; xhat];
        solution_found = i >= iterations;
        
        % Updates
        fis = [fis; oracle(xhat)];
        m_1 = m;
        deltax_1 = xhat - xhat_1;
        xhat_1 = xhat;    
    end
end