function [x_star, xhats, xs, vs, ms, Es, deltas, zbars] = ...
            drecexpbary_2nd_custom(oracle, m0, xhat0, vhat0, nu, sigma, ...
                                   lambda, gamma_, damp_fac, tau, ...
                                   iterations, options)
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

    MAX_ITERATIONS = 100;
    
    dt = 1e-3;
    gamma = 1;
    b = 0.01;
    
    % Previous values
    m_1 = m0;
    x_1 = xhat0;
    xhat_1 = xhat0;
    v_1 = vhat0;
    E_1 = m0*(vhat0'*vhat0 + gamma_*oracle(xhat0));
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    delta_vhat_1 = zeros(size(vhat0));
        
    % Accumulated values
    xhats = xhat0';
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
       %zbar = accel_fun(m_1, vhat_1, delta_vhat_1);
       
       zbar = zeros(size(vhat0));
       zbar = double(zbar);
       sigma = double(sigma);
       
       z = zeros(length(zbar), 1);
       for j = 1:length(z)
           z(j) = gaussianrnd(zbar(j), sigma);
       end
    
       x_n = xhat_1 + tau*v_1;
       e_n = exp(-nu*oracle(x_n));
       m_n = lambda*m_1 + e_n;
       xhat_n = (1/m_n)*(lambda*m_1*xhat_1 + e_n*x_n);
       norm_v_n = sqrt((1/m_n)*(m_1*v_1.'*v_1 + ...
                               gamma_*m_1*oracle(x_1) + ...
                               - damp_fac*m_1*v_1.'*v_1 - ...
                               gamma_*m_n*oracle(x_n)));
       norm_v_1_z = sqrt((v_1 + z).'*(v_1 + z));
       v_n = (v_1 + z)*(norm_v_n/norm_v_1_z);
       
       % Current value of position
       solution_found = i >= iterations;
              
       % Updates
       m_1 = m_n;
       x_1 = x_n;
       delta_xhat_1 = xhat_n - xhat_1;
       xhat_1 = xhat_n;
       v_1 = v_n;
       
       xs = [xs; x_n'];
       vs = [vs; v_n'];
       xhats = [xhats; xhat_n'];
       ms = [ms; m_n];
       deltas = [deltas; delta_xhat_1'];
       Es = [Es; oracle(x_n)];
       
       i = i + 1;
       
       if(options.verbose)
          coord_string_xhat = vec2str(xhat_n);
          coord_string_v = vec2str(v_n);
    
          args = [i, xhat_n', v_n', oracle(x_n)];
          disp(sprintf([' iter = %3d: ', ...
                        ' xhat: ', '(', coord_string_xhat, ')', ...
                        ' v: ', '(', coord_string_v, ')', ...
                        ' fval: %.3f'], args));
       end
       
       wb = wb.update_waitbar(i, MAX_ITERATIONS);
    end
    
    x_star = xhats(end, :);
end
