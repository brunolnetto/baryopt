function [x, xs] = recexpbary(oracle, x0, nu, ...
                              sigma, zeta, lambda, ...
                              iterations)
    xhat_1 = x0;
    m_1 = 0;
    deltax_1 = zeros(size(x0));
    solution_found = false;    
    
    xs = [];
    
    i = 1;
    while(~solution_found)
        z = normrnd(zeta*deltax_1, sigma);
        x = xhat_1 + z;
        
        ei = exp(-nu*oracle(x));
        
        m = lambda*m_1 + ei;
        xhat = (1/m)*(lambda*m_1*xhat_1 + x*ei);
        
        xs = [xs; xhat];
        norm(xhat - xhat_1)
        solution_found = i >= iterations;
        
        % Updates
        i = i + 1;
        m_1 = m;
        deltax_1 = xhat - xhat_1;
        xhat_1 = xhat;    
    end
end