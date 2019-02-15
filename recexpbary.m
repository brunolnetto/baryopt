function [x, xs] = recexpbary(oracle, x0, nu, ...
                              sigma, zeta, lambda, gamma, ...
                              iterations, method)
% Recursive barycenter algorithm for direct optimization
% https://arxiv.org/abs/1801.10533
% In:
%   - oracle [function]: Oracle function
%   - x0 []: Query values
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
        
        if(strcmp(method, 'shape'))
%             an = rand(size(x0));
            an = normrnd(zeta*deltax_1, sigma);
            z = an*(fis(i-1)/max(fis))^gamma;
        elseif(strcmp(method, 'normal'))
            z = normrnd(zeta*deltax_1, sigma);
        else
            error('method MUST be rather normal or shape');
        end
        x = xhat_1 + z;
        
        ei = exp(-nu*oracle(x));
        
        m = lambda*m_1 + ei;
        xhat = (1/m)*(lambda*m_1*xhat_1 + x*ei);
        
        xs = [xs; xhat];
        solution_found = i >= iterations;
        
        % Updates
        fis = [fis; oracle(xhat)];
        m_1 = m;
        deltax_1 = xhat - xhat_1;
        xhat_1 = xhat;    
    end
end