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

x0 = init_val*ones(n, 1);
m0 = exp(-nu*oracle(x0));

iterations = 20;
n_iterations = 1000;

wb = my_waitbar('Calculating minimum...');
lambda_zs = [1 0.66 0.33];

zbars_total = zeros(iterations, 2);

% Reference values
x_test = {};
for i = 1:n_iterations
    [x, xs] = drecexpbary_custom(oracle, m0, x0, ...
                                 nu, sigma, lambda, 0, ...
                                 iterations, false);

    x_test{end+1} = xs;

    wb.update_waitbar(i, n_iterations);
end

x_acc = zeros(size(x_test{1}));
for i = 1:length(x_test)
    x_i = x_test{i};
    x_acc = x_acc + x_test{i};
end    

xs_lambdaz0 = x_acc/n_iterations;

% Accelerated values
zbars_lambdaz = {};
xs_lambdaz = {};
for lambda_z = lambda_zs
    zbars_total = zeros(iterations, 2);
    x_test = {};
    
    for i = 1:n_iterations
        [x, xs] = drecexpbary_custom(oracle, m0, x0, ...
                                     nu, sigma, lambda, lambda_z, ...
                                     iterations, true);
        
        x_test{end+1} = xs;
        zbars_total = zbars_total + zbars;

        wb.update_waitbar(i, n_iterations);
    end
    
    x_acc = zeros(size(x_test{1}));
    for i = 1:length(x_test)
        x_i = x_test{i};
        x_acc = x_acc + x_test{i};
    end    
    
    xs_lambdaz{end+1} = x_acc/n_iterations;    
    zbars_lambdaz{end+1} = zbars_total/iterations;
end

delta_xlambdaz = {};
for i = 1:length(xs_lambdaz)
    xs_lamb = xs_lambdaz{i};
    
    delta_xjs = [];
    for j = 2:length(xs_lamb)
        delta_xj = xs_lamb(j, :) - xs_lamb(j-1, :);
        delta_xjs = [delta_xjs; delta_xj];
    end
    
    delta_xlambdaz{end+1} = delta_xjs;
end

delta_xlambdaz0 = [];
for j = 2:length(xs_lamb)
    delta_xj = xs_lambdaz0(j, :) - xs_lambdaz0(j-1, :);
    delta_xlambdaz0 = [delta_xlambdaz0; delta_xj];
end

% Delta x plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\frac{\Delta \hat{x}_1(\lambda_z) - \Delta \hat{x}_1(0)}{\Delta \hat{x}_1(0)}$', ...
                       '$\frac{\Delta \hat{x}_2(\lambda_z) - \Delta \hat{x}_1(0)}{\Delta \hat{x}_1(0)}$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'0.33', '0.66', '1'}, ...
                       {'0.33', '0.66', '1'}};
plot_config.pos_multiplots = [1, 1, 2, 2];
plot_config.markers = {{'k-', 'r--', '.b-'}, ...
                       {'k-', 'r--', '.b-'}};

delta_first = [(delta_xlambdaz{3}(:, 1) - delta_xlambdaz0(:, 1))./delta_xlambdaz0(:, 1), ...
               (delta_xlambdaz{3}(:, 2) - delta_xlambdaz0(:, 2))./delta_xlambdaz0(:, 2)];
delta_multi = [(delta_xlambdaz{2}(:, 1) - delta_xlambdaz0(:, 1))./delta_xlambdaz0(:, 1), ...
               (delta_xlambdaz{1}(:, 1) - delta_xlambdaz0(:, 1))./delta_xlambdaz0(:, 1), ...
               (delta_xlambdaz{2}(:, 2) - delta_xlambdaz0(:, 2))./delta_xlambdaz0(:, 2), ...
               (delta_xlambdaz{1}(:, 2) - delta_xlambdaz0(:, 2))./delta_xlambdaz0(:, 2)];

iters = 1:length(delta_first(:, 1));
hfigs_delta = my_plot(iters, {delta_first, delta_multi}, plot_config);

% Acceleration plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'0', '0.33', '0.66', '1'}, ...
                       {'0', '0.33', '0.66', '1'}};
plot_config.pos_multiplots = [1, 1, 1, 2, 2, 2];
plot_config.markers = {{'k-', 'r--', '.b-', 'sg'}, ...
                       {'k-', 'r--', '.b-', 'sg'}};

n = length(zbars_lambdaz{3}(:, 1));
zbars_first = [zeros(n, 1), zeros(n, 1)];
zbars_multi = [zbars_lambdaz{3}(:, 1), ...
               zbars_lambdaz{2}(:, 1), ...
               zbars_lambdaz{1}(:, 1), ...
               zbars_lambdaz{3}(:, 2), ...
               zbars_lambdaz{2}(:, 2), ...
               zbars_lambdaz{1}(:, 2)];

iters = 1:length(zbars(:, 1));
hfigs_z = my_plot(iters, {zbars_first, zbars_multi}, plot_config);

% x-coordinate plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\frac{\hat{x}_1(\lambda_z) - \hat{x}_1(0)}{\hat{x}_1(0)}$', ...
                       '$\frac{\hat{x}_2(\lambda_z) - \hat{x}_2(0)}{\hat{x}_2(0)}$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'0.33', '0.66', '1'}, {'0.33', '0.66', '1'}};
plot_config.pos_multiplots = [1, 1, 2, 2];
plot_config.markers = {{'r--', '.b-', 'sg'}, {'r--', '.b-', 'sg'}};

x_first = [(xs_lambdaz{3}(:, 1) - xs_lambdaz0(:, 1))./(xs_lambdaz0(:, 1)), ...
           (xs_lambdaz{3}(:, 2) - xs_lambdaz0(:, 2))./xs_lambdaz0(:, 2)];
x_multi = [(xs_lambdaz{2}(:, 1) - xs_lambdaz0(:, 1))./xs_lambdaz0(:, 1), ...
           (xs_lambdaz{1}(:, 1) - xs_lambdaz0(:, 1))./xs_lambdaz0(:, 1), ...
           (xs_lambdaz{2}(:, 2) - xs_lambdaz0(:, 2))./xs_lambdaz0(:, 2), ...
           (xs_lambdaz{1}(:, 2) - xs_lambdaz0(:, 2))./xs_lambdaz0(:, 2)];

iters = 1:length(xs_lambdaz{1}(:, 1));
[hfigs_x_x0, axs] = my_plot(iters, {x_first, x_multi}, plot_config);

axs{1}{1}.XLim = [0 20];
axs{1}{2}.XLim = [0 20];

% Save folder
path = [pwd '/../imgs/'];
fname = ['xbars_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_delta, [path, fname], 'epsc');

path = [pwd '/../imgs/'];
fname = ['xbars_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_z, [path, fname], 'epsc');

path = [pwd '/../imgs/'];
fname = ['xbars_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_x_x0, [path, fname], 'epsc');


