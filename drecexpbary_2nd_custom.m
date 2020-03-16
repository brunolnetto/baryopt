function [x_star, xhats, ...
          xs, vhats, vs, ms, ...
          Es, deltas, zbars] = ...
            drecexpbary_2nd_custom(oracle, m0, xhat0, vhat0, nu, ...
                                   sigma, lambda, iterations, accel_fun)
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
    MAX_ITERATIONS = 100;
    
    dt = 1e-3;
    gamma = 1;
    b = 0.01;
    
    % Previous values
    F = 0;
    m_1 = m0;
    xhat_1 = xhat0;
    vhat_1 = vhat0;
    v_1 = vhat0;
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    delta_vhat_1 = zeros(size(vhat0));
        
    % Accumulated values
    xhats = xhat0';
    vhats = vhat0';
    zbars = [];
    ms = m0;
    deltas = [];
    xs = [];
    vs = [];
    Es = [];
    
    solution_found = false;
    
    wb = my_waitbar('Calculating minimum...');
    
    i = 1;
    while(~solution_found)
       % Calculation of mean for stochastic function
       zbar = accel_fun(m_1, vhat_1, delta_vhat_1);
       
       zbar = double(zbar);
       sigma = double(sigma);
       
       z = zeros(length(zbar), 1);
       for j = 1:length(z)
           z(j) = gaussianrnd(zbar(j), sigma);
       end
       
       % Current value of position
       v = vhat_1 + z;
       x = xhat_1 + dt*v_1;
       
       E_i = m_1*(oracle(x) + gamma*v'*v);
       
       % Mass component
       e_i = exp(-nu*E_i);
       m = lambda*m_1 + e_i;
       
       % Barycenter point
       vhat = (1/m)*(m_1*vhat_1 + e_i*v);
       xhat = (1/m)*(m_1*xhat_1 + e_i*x);
       
       solution_found = i >= iterations;
              
       % Updates
       m_1 = m;
       delta_xhat_1 = xhat - xhat_1;
       xhat_1 = xhat;
       vhat_1 = vhat;
       v_1 = v;
       
       xs = [xs; x'];
       vs = [vs; v'];
       
       xhats = [xhats; xhat'];
       vhats = [vhats; vhat'];
       
       ms = [ms; m];
       
       zbars = [zbars; zbar'];
       deltas = [deltas; delta_xhat_1'];
       
       Es = [Es; E_i];
       
       i = i + 1;
       
       wb = wb.update_waitbar(i, MAX_ITERATIONS);
       
       if(i == MAX_ITERATIONS)
           MAX_ITERATIONS = MAX_ITERATIONS*2;
       end
    end
    
    x_star = xhats(end, :);
end
