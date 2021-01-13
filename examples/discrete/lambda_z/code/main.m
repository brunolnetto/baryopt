close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 1;
lambda = 1;

% Recursive version
syms x y

x_0 = 1;
y_0 = -1;

func = (x - x_0)^2 + (y - y_0)^2;
oracle = @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;

x0 = init_val*zeros(n, 1);
m0 = 0;

iterations = 25;
n_iterations = 150;

wb = my_waitbar('Calculating minimum...');
lambda_zs = {1, 0.9 0.8, 0.7};

% Accelerated values
zbars_lambdaz = {};
xs_lambdaz = {};
deltas_lambdaz = {};

idx_1 = 1;
for lambda_z = lambda_zs
    lambda_z = lambda_z{1};
    
    zbars_acc = zeros(iterations, 2);
    x_acc = zeros(iterations, 2);
    deltas_acc = zeros(iterations-1, 2);
    
    for i = 1:n_iterations
        accel_fun = @(m_1, xhat_1, delta_xhat_1) ...
            integrated_accel(m_1, xhat_1, delta_xhat_1, ...
                             lambda_z, nu, oracle);
        
        [~, xhats, ~, ~, deltas, zbars, ~, ~] = ...
            drecexpbary_custom(oracle, m0, x0, ...
                               nu, sigma, ...
                               lambda, iterations, ...
                               accel_fun, ...
                               struct('verbose', true));
        
        clear(func2str(@integrated_accel));
        
        xhats = xhats(1:end-1, :);
        zbars = zbars(1:end-1, :);
        deltas = deltas(1:end-1, :);
        
        zbars_acc = zbars_acc + zbars;
        x_acc = x_acc + xhats;
        deltas_acc = deltas_acc + deltas;

        wb.update_waitbar(i + (idx_1-1)*n_iterations, ...
                          length(lambda_zs)*n_iterations);
    end
        
    xs_lambdaz{end+1} = x_acc/n_iterations;
    zbars_lambdaz{end+1} = zbars_acc/iterations;
    deltas_lambdaz{end+1} = deltas_acc/iterations;
    idx_1 = idx_1 + 1;
end

lambda_zs = cellfun(@num2str, lambda_zs, 'UniformOutput', false);
markers_ = {'k-', 'r--', '.b-', 'sg'};
markers = {markers_{1:length(lambda_zs)}};

% Delta x plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\Delta \hat{x}_1(\lambda_z)$', ...
                       '$\Delta \hat{x}_2(\lambda_z)$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {lambda_zs, lambda_zs};
plot_config.pos_multiplots = [repmat(1, 1, length(lambda_zs)-1), ...
                              repmat(2, 1, length(lambda_zs)-1)];
plot_config.markers = {markers, markers};
plot_config.axis_style = 'square';

delta_first = [deltas_lambdaz{1}(:, 1), deltas_lambdaz{1}(:, 2)];
delta_multi = [deltas_lambdaz{2}(:, 1), deltas_lambdaz{3}(:, 1), ...
               deltas_lambdaz{4}(:, 1), deltas_lambdaz{2}(:, 2), ...
               deltas_lambdaz{3}(:, 2), deltas_lambdaz{4}(:, 2)];

deltas = {delta_first, delta_multi};

iters = 1:length(delta_first(:, 1));
hfigs_delta = my_plot(iters, deltas, plot_config);

% Acceleration plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {lambda_zs, lambda_zs};
plot_config.pos_multiplots = [repmat(1, 1, length(lambda_zs)-1), ...
                              repmat(2, 1, length(lambda_zs)-1)];
plot_config.markers = {markers, markers};
plot_config.axis_style = 'square';

n = length(zbars_lambdaz{1}(:, 1));
zbars_first = [zbars_lambdaz{1}(:, 1), zbars_lambdaz{1}(:, 2)];
zbars_multi = [zbars_lambdaz{2}(:, 1), zbars_lambdaz{3}(:, 1), ...
               zbars_lambdaz{4}(:, 1), zbars_lambdaz{2}(:, 2), ...
               zbars_lambdaz{3}(:, 2), zbars_lambdaz{4}(:, 2)];

zs = {zbars_first, zbars_multi};
           
iters = 1:length(zbars(:, 1));
hfigs_z = my_plot(iters, zs, plot_config);

% x-coordinate plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\hat{x}_1(\lambda_z)$', ...
                       '$\hat{x}_2(\lambda_z)$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {lambda_zs, lambda_zs};
plot_config.pos_multiplots = [repmat(1, 1, length(lambda_zs)-1), ...
                              repmat(2, 1, length(lambda_zs)-1)];
plot_config.markers = {markers, markers};
plot_config.axis_style = 'square';

x_first = [xs_lambdaz{1}(:, 1), xs_lambdaz{1}(:, 2)];
x_multi = [xs_lambdaz{2}(:, 1), xs_lambdaz{3}(:, 1), ...
           xs_lambdaz{4}(:, 1), xs_lambdaz{2}(:, 2), ...
           xs_lambdaz{3}(:, 2), xs_lambdaz{4}(:, 2)];

xs = {x_first, x_multi};

iters = 1:length(xs_lambdaz{1}(:, 1));
[hfigs_x_x0, axs] = my_plot(iters, xs, plot_config);

axs{1}{1}.XLim = [0 iterations];
axs{1}{2}.XLim = [0 iterations];

% Save folder
path = [pwd '/../imgs/'];
fname = ['delta_', ...
         sprintf('lamb%.2f', 100*lambda), ...
         sprintf('sigma%.2f', 100*sigma), ...
         sprintf('nu%d', nu), ...
         sprintf('lambz%.2f', 100*lambda_z)];
saveas(hfigs_delta, [path, fname], 'epsc');

fname = ['zbar_', ...
         sprintf('lamb%.2f', 100*lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_z, [path, fname], 'epsc');

fname = ['xbars_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_x_x0, [path, fname], 'epsc');


