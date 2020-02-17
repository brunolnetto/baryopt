close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 0.5;
lambda = 1;

% Recursive version
x_0 = 1;
y_0 = -1;
oracle = @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;

syms x y
func = (x - 1)^2 + (y + 1)^2;
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
m0 = exp(-nu*oracle(x0));

iterations = 20;
n_iterations = 2000;

x_test = {};

wb = my_waitbar('Calculating minimum...');
lambda_zs = [1 0.66 0.33];

zbars_lambdaz = {};

for lambda_z = lambda_zs
    zbars_total = zeros(iterations, 2);
    for i = 1:n_iterations
        [x, xs] = drecexpbary_custom(oracle, m0, x0, ...
                                     nu, sigma, lambda, lambda_z, ...
                                     iterations, true);
        x_test{end+1} = xs;

        zbars_total = zbars_total + zbars;

        wb.update_waitbar(i, n_iterations);
    end

    zbars_lambdaz{end+1} = zbars_total/iterations;
end

% Acceleration plot
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'1', '0.66', '0.33'}, {'1', '0.66', '0.33'}};
plot_config.pos_multiplots = [1, 1, 2, 2];
plot_config.markers = {{'k-', 'r--', '.b-'}, {'k-', 'r--', 'b.-'}};

zbars_first = [zbars_lambdaz{1}(:, 1), zbars_lambdaz{1}(:, 2)];
zbars_multi = [zbars_lambdaz{2}(:, 1), zbars_lambdaz{3}(:, 1), ...
               zbars_lambdaz{2}(:, 2), zbars_lambdaz{3}(:, 2)];

iters = 1:length(zbars(:, 1));
hfigs_s = my_plot(iters, {zbars_first, zbars_multi}, plot_config);

% Save folder
path = [pwd '/../imgs/'];
fname = ['zbars_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_s, [path, fname], 'epsc');

