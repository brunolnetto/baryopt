close all
clear all
clc

n = 2;
init_val = 0;

% Method hyperparameter
nu = 10;
sigma = 1;
Sigma = sigma*eye(n);
lambda = 1;
lambda_z = 1;
zeta = 2;

% Recursive version
syms x y
x_0 = 1;
y_0 = -1;
func = (x - x_0)^2 + (y - y_0)^2;
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

iterations = 100;
n_iterations = 100;

wb = my_waitbar('Calculating minimum...');
is_accel = false;

accel_fun_nonaccel = @(m_1, xhat_1, delta_xhat_1) ...
                     non_accel(m_1, xhat_1, delta_xhat_1);
accel_fun_pait = @(m_1, xhat_1, delta_xhat_1) ...
                  pait_accel(m_1, xhat_1, delta_xhat_1, zeta);
accel_fun_integrated = @(m_1, xhat_1, delta_xhat_1) ...
                            integrated_accel(m_1, xhat_1, ...
                                                  delta_xhat_1, lambda_z, ...
                                                  nu, oracle);

% Methods and names
accel_funs = {accel_fun_nonaccel, ...
              accel_fun_pait, ...
              accel_fun_integrated};
names = {'Not accelerated', ...
         'Pait acceleration', ...
         'Integrated acceleration'};

% Mean value for methods
xhat_means = {};
deltas_means = {};
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
    ms_tests = {};
    
    % Method understood
    for j = 1:n_iterations    
        [x_star, xhats, xs, m, deltas, zbars, Fs_n, Fbars_n] = ...
         drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
                            lambda, iterations, accel_fun);

        xstar_tests{end+1} = x_star;
        xs_tests{end+1} = xs;
        xhat_tests{end+1} = xhats;
        delta_tests{end+1} = deltas;
        zbars_tests{end+1} = zbars;
        ms_tests{end+1} = m;

        wb.update_waitbar(j + (i-1)*n_iterations, length(names)*n_iterations);
    end
    
    % Take the mean value
    
    % X star
    xstar = zeros(1, 2);
    for i = 1:length(xhat_tests)
        x_aux = xstar_tests{i};
        xstar = x_star + x_aux;
    end

    xstar_mean = xstar/n_iterations;
    
    % X hat and x
    xhat_mean = zeros(size(xhat_tests{1}));
    x_mean = zeros(size(xs_tests{1}));
    for i = 1:length(xhat_tests)
        x_mean = x_mean + xs_tests{i};
        xhat_mean = xhat_mean + xhat_tests{i};
    end

    xhat_mean = xhat_mean/n_iterations;
    x_mean = x_mean/n_iterations;
    
    % m
    m_mean = zeros(size(ms_tests{1}));
    for i = 1:length(xhat_tests)
        m_mean = m_mean + ms_tests{i};
        m_mean = m_mean + ms_tests{i};
    end
    
    m_mean = m_mean/n_iterations;
    
    % Deltas
    delta_mean = zeros(size(delta_tests{1}));
    for i = 1:length(zbars_tests)
        delta_mean = delta_mean + delta_tests{i};
    end

    delta_mean = delta_mean/n_iterations;
    
    % zbars
    zbar_mean = zeros(size(zbars_tests{1}));
    for i = 1:length(zbars_tests)
        zbar_mean = zbar_mean + zbars_tests{i};
    end

    zbar_mean = zbar_mean/length(zbars_tests);
    
    % function value
    fs_mean = [];
    for i = 1:iterations
        f_i = oracle(xhat_mean(i, :));
        fs_mean = [fs_mean; f_i];
    end
    
    xstar_means{end+1} = xstar_mean;
    x_means{end+1} = x_mean;
    xhat_means{end+1} = xhat_mean;
    zbars_means{end+1} = zbar_mean;
    ms_means{end+1} = m_mean;
    fs_means{end+1} = fs_mean;
    deltas_means{end+1} = delta_mean;
end

% -------- Plots --------
title = sprintf('$\\nu = %.2f$, $\\sigma = %.2f$, $\\zeta = %.2f$, $\\lambda_z = %.2f$', ...
                nu, sigma, zeta, lambda_z);
legend_atom = {'$\bar{z} = 0$', ...
               '$\bar{z} = \zeta \, \Delta \hat{x}$', ...
               '$\bar{z} = barycenter([\Delta \hat{x}_1, \cdots, \Delta \hat{x}_{n-1}])$'};

markers_atom = {'-', '--', '.-'};

% zbar plot
plot_config_z.titles = {'', ''};
plot_config_z.xlabels = {'', 'Iterations'};
plot_config_z.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config_z.grid_size = [2, 1];
plot_config_z.plot_type = 'stem';
plot_config_z.legends = {legend_atom, legend_atom};
plot_config_z.pos_multiplots = [1, 1, 2, 2];
plot_config_z.markers = {markers_atom, markers_atom};

zbars_means23 = double([zbars_means{2}, zbars_means{3}]);
zbars_ = {double(zbars_means{1}), ...
          [zbars_means23(:, 1), zbars_means23(:, 3), ...
           zbars_means23(:, 2), zbars_means23(:, 4)]};

[m, n] = size(zbars(:, 1));
iters = (1:length(zbars(:, 1)))';
hfigs_zmean = my_plot(iters, zbars_, plot_config_z);

% function plot
fs_ = {double(fs_means{1}), ...
       double([fs_means{2}, fs_means{3}])};

plot_config_f.titles = {[title, ' $f(x, y) := ', latex(func),'$']};
plot_config_f.xlabels = {'Iterations'};
plot_config_f.ylabels = {'$f(x, y)$'};
plot_config_f.grid_size = [1, 1];
plot_config_f.legends = legend_atom;
plot_config_f.pos_multiplots = [1, 1];
plot_config_f.markers = markers_atom;

[m, n] = size(fs_{1});
iters = 1:m;
hfigs_fs = my_plot(iters, fs_, plot_config_f);

% Delta xhat plot
deltas_means23 = double([deltas_means{2}, deltas_means{3}]);
delta_xhat_ = {double(deltas_means{1}), ...
               double([deltas_means23(:, 1), deltas_means23(:, 3), ...
                       deltas_means23(:, 2), deltas_means23(:, 4)])};

plot_config_delta.titles = {title, ''};
plot_config_delta.xlabels = {'', 'Iterations'};
plot_config_delta.ylabels = {'$\Delta \hat{x}_1$', '$\Delta \hat{x}_2$'};
plot_config_delta.grid_size = [2, 1];
plot_config_delta.plot_type = 'stem';
plot_config_delta.legends = {legend_atom, legend_atom};
plot_config_delta.pos_multiplots = [1, 1, 2, 2];
plot_config_delta.markers = {markers_atom, markers_atom};

iters = (1:length(delta_xhat_{1}))';
hfigs_delta_mean = my_plot(iters, delta_xhat_, plot_config_delta);

% xhat plot
xhat_means23 = double([xhat_means{2}, xhat_means{3}]);
xhat_ = {double(xhat_means{1}), ...
          double([xhat_means23(:, 1), xhat_means23(:, 3), ...
                  xhat_means23(:, 2), xhat_means23(:, 4)])};''

plot_config_xhat.titles = {title, ''};
plot_config_xhat.xlabels = {'', 'Iterations'};
plot_config_xhat.ylabels = {'$\hat{x}_1$', '$\hat{x}_2$'};
plot_config_xhat.grid_size = [2, 1];
plot_config_xhat.plot_type = 'stem';
plot_config_xhat.legends = {legend_atom, legend_atom};
plot_config_xhat.pos_multiplots = [1, 1, 2, 2];
plot_config_xhat.markers = {markers_atom, markers_atom};

iters = (1:length(xhat_{1}))';
hfigs_xmean = my_plot(iters, xhat_, plot_config_xhat);

% x plot
x_means23 = double([x_means{2}, x_means{3}]);
x_ = {double(x_means{1}), double([x_means23(:, 1), x_means23(:, 3), ...
                          x_means23(:, 2), x_means23(:, 4)])};

plot_config_xhat.titles = {title, ''};
plot_config_xhat.xlabels = {'', 'Iterations'};
plot_config_xhat.ylabels = {'$x_1$', '$x_2$'};
plot_config_xhat.grid_size = [2, 1];
plot_config_xhat.plot_type = 'stem';
plot_config_xhat.legends = {legend_atom, legend_atom};
plot_config_xhat.pos_multiplots = [1, 1, 2, 2];
plot_config_xhat.markers = {markers_atom, markers_atom};

iters = (1:length(xhat_{1}))';
hfigs_xmean = my_plot(iters, xhat_, plot_config_xhat);

% m plot
m_ = {double(ms_means{1}), double([ms_means{2}, ms_means{3}])};

plot_config_m.titles = {''};
plot_config_m.xlabels = {'Iterations'};
plot_config_m.ylabels = {'$x_1$'};
plot_config_m.grid_size = [1, 1];
plot_config_m.plot_type = 'stem';
plot_config_m.legends = {legend_atom};
plot_config_m.pos_multiplots = [1, 1];
plot_config_m.markers = {markers_atom};

iters = (1:length(m_{1}))';
hfigs_xmean = my_plot(iters, m_, plot_config_m);

