close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 10;
sigma = 0.5;
Sigma = sigma*eye(n);
lambda = 1;
lambda_z = 0.3;
zeta = 0.5;

% Recursive version
syms x y
x_0 = 1;
y_0 = -1;
func = (x - x_0)^2 + (y - y_0)^2;
grad_f = gradient(func);
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

% syms x y
% func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
%               - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
%               - 1/3*exp(-(x+1)^2 - y^2);
% 
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = 1;

iterations = 10;
n_iterations = 100;

wb = my_waitbar('Calculating minimum...');
is_accel = false;

accel_fun_nonaccel = @(m_1, xhat_1, delta_xhat_1) ...
                     non_accel(m_1, xhat_1, delta_xhat_1);
accel_fun_pait = @(m_1, xhat_1, delta_xhat_1) ...
                  pait_accel(m_1, xhat_1, delta_xhat_1, zeta);
accel_fun_integrated = @(m_1, xhat_1, delta_xhat_1) ...
                      integrated_accel(m_1, xhat_1, delta_xhat_1, ...
                                       lambda_z, nu, oracle);

% accel_funs = {accel_fun_nonaccel, accel_fun_pait, accel_fun_integrated};
% names = {'Not accelerated', ...
%          'Pait acceleration', ...
%          'Integrated acceleration'};

accel_funs = {accel_fun_nonaccel, accel_fun_pait, accel_fun_integrated};
names = {'Not accelerated', ...
         'Pait acceleration', ...
         'Integrated acceleration'};

xhat_means = {};
x_means = {};
zbars_means = {};
xstar_means = {};
ms_means = {};
fs_means = {};
for i = 1:length(names)
    accel_fun = accel_funs{i};
    name = names{i};
    
    xstar_tests = {};
    xs_tests = {};
    xhat_tests = {};
    delta_tests = {};
    zbars_tests = {}; 
    F_tests = {};
    ms_tests = {};
    Fbar_tests = {}; 
    for i = 1:n_iterations    
        [x_star, xhats, xs, m, deltas, zbars, Fs_n, Fbars_n] = ...
         drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
                            lambda, iterations, accel_fun);

        xstar_tests{end+1} = x_star;
        xs_tests{end+1} = xs;
        xhat_tests{end+1} = xhats;
        delta_tests{end+1} = deltas;
        zbars_tests{end+1} = zbars;
        ms_tests{end+1} = m;

        wb.update_waitbar(i, n_iterations);
    end

    xstar = zeros(1, 2);
    for i = 1:length(xhat_tests)
        x_aux = xstar_tests{i};
        xstar = x_star + x_aux;
    end

    xstar_mean = xstar/n_iterations;

    xhat_mean = zeros(size(xhat_tests{1}));
    x_mean = zeros(size(xs_tests{1}));
    for i = 1:length(xhat_tests)
        x_mean = x_mean + xs_tests{i};
        xhat_mean = xhat_mean + xhat_tests{i};
    end

    xhat_mean = xhat_mean/n_iterations;
    x_mean = x_mean/n_iterations;

    m_mean = zeros(size(ms_tests{1}));
    for i = 1:length(xhat_tests)
        m_mean = m_mean + ms_tests{i};
        m_mean = m_mean + ms_tests{i};
    end
    
    m_mean = m_mean/n_iterations;
    
    zbar_mean = zeros(size(zbars_tests{1}));
    for i = 1:length(zbars_tests)
        zbar_mean = zbar_mean + zbars_tests{i};
    end

    zbar_mean = zbar_mean/length(zbars_tests);

    fs_mean = [];
    for i = 1:iterations
        f_i = oracle(x_mean(i, :));
        fs_mean = [fs_mean; f_i];
    end
    
    xstar_means{end+1} = xstar_mean;
    x_means{end+1} = x_mean;
    xhat_means{end+1} = xhat_mean;
    zbars_means{end+1} = zbar_mean;
    ms_means{end+1} = m_mean;
    fs_means{end+1} = fs_mean;
end

% plot_config.titles = {'', ''};
% plot_config.xlabels = {'', 'Iterations'};
% plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
% plot_config.grid_size = [2, 1];
% plot_config.plot_type = 'stem';
% 
% iters = 1:length(zbars(:, 1));
% hfigs_zmean = my_plot(iters, zbar_mean, plot_config);
% 
% plot_config.titles = {['$f(x, y) := ', latex(func),'$']};
% plot_config.xlabels = {'Iterations'};
% plot_config.ylabels = {'$f(x, y)$'};
% plot_config.grid_size = [1, 1];
% plot_config.plot_type = 'stem';
% 
% iters = 1:length(fs);
% hfigs_fs = my_plot(iters, fs, plot_config);
% 
% plot_config.titles = {'', ''};
% plot_config.xlabels = {'', 'Iterations'};
% plot_config.ylabels = {'Amplitude', 'Amplitude'};
% plot_config.grid_size = [2, 1];
% plot_config.plot_type = 'stem';
% plot_config.legends = {{'$\hat{x}_1$', '$x$'}, ...
%                        {'$\hat{x}_2$', '$x$'}};
% plot_config.pos_multiplots = [1, 2];
% plot_config.markers = {{'-', '--'}, {'-', '--'}};
% 
% ys = {xhat_mean, x_mean};
% 
% iters = 1:length(xs(:, 1));
% hfigs_xmean = my_plot(iters, xhat_mean, plot_config);
% 
