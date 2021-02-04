function [x_star, xhats, xs, ms, ...
          deltas, zbars, Fs_n, Fbars_n] = ...
            drecexpbary_custom(oracle, m0, xhat0, nu, sigma, ...
                               lambda, iterations, accel_fun, options)
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
    
    if(nargin == 8)
        options = struct('verbose', false, ...
                         'is_progress_visible', false);
    end
    
    if(~isfield(options, 'verbose'))
        options.verbose = false;
    end
    
    if(~isfield(options, 'is_progress_visible'))
        options.is_progress_visible = false;
    end
    
    alpha_n = 1;
    deltas = [];
        
    % Mass components
    m_1 = terop(m0 ~= 0, m0, eps);
    xhat_1 = xhat0;
    x_1 = xhat0;
    zbar0 = zeros(size(xhat0));
    z_bar_1 = zeros(size(xhat0));
    delta_xhat_1 = zeros(size(xhat0));
    
    xhats = xhat0';
    ms = m0;
    deltas = [];
    zbars = zeros(size(xhat0'));
    Fs_n = [];
    Fbars_n = [];
    xs = [];
    fs = oracle(xhat0);
    
    solution_found = false;
    
    coord_string_xhat = vec2str(xhat0);
    coord_string_zbar = vec2str(zbar0);
    
    xhat = xhat0;
    zbar = z_bar_1;
    max_grad0 = 0;
    
    j = 1;
    
    if(options.is_progress_visible)
        wb = my_waitbar('Calculating minimum...');
    end
    
    while(~solution_found)
       if(options.verbose)
           args = [j, xhat', zbar', oracle(xhat)];
           disp(sprintf([' iter = %3d: ', ...
                         ' xhat: ', '(', coord_string_xhat, ')', ...
                         ' zbar: ', '(', coord_string_zbar, ')', ...
                         ' fval: %.3f'], args));
       end
        
       % Expected value for stochastic function
       zbar = accel_fun(m_1, x_1, delta_xhat_1);
       
       z = zeros(length(zbar), 1);
       for k = 1:length(zbar)            
           sigma_n = double(sigma);
           
           z(k) = gaussianrnd(double(zbar(k)), sigma_n);
       end
       
       % Current value of position
       x = xhat_1 + z;
       
       % Mass component
       fx = oracle(x);
       
       e_i = exp(-nu*fx);
       m = lambda*m_1 + e_i;
       
       F_n = vpa(e_i/m);
       Fbar_n = (m_1/m)*F_n;
       
       % Barycenter point
       xhat = vpa(xhat_1 + F_n*z);              
       
       solution_found = j >= iterations;
       
       % Updates
       m_1 = vpa(m);
       delta_xhat_1 = vpa(xhat - xhat_1);
       xhat_1 = vpa(xhat);
       x_1 = vpa(x);
       
       fs = [fs; fx];
       xs = [xs; x'];
       xhats = [xhats; xhat'];
       ms = [ms; m];
       deltas = [deltas; delta_xhat_1'];
       Fs_n = [Fs_n; F_n];
       Fbars_n = [Fbars_n; Fbar_n];
       zbars = [zbars; zbar'];
       
       if(options.is_progress_visible)
           wb = wb.update_waitbar(j, iterations);
       end
       
       j = j + 1;
    end
    
    if(options.is_progress_visible)
       wb.close_window();
   end
    
    
    x_star = xhats(end, :)';
end
