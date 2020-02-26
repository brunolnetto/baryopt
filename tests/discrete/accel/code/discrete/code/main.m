close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 0.5;
lambda = 1;
lambda_z = 0.1;

% Recursive version
syms x y
x_0 = 1;
y_0 = -1;
oracle = @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;
func = (x - x_0)^2 + (y - y_0)^2;

x0 = init_val*ones(n, 1);
m0 = 0;

iterations = 50;
is_accel = true;

x_tests = {};
[x, xs, m, deltas, zbars]...
    = drecexpbary_custom(oracle, m0, x0, nu, sigma, lambda, ...
                         lambda_z, iterations, is_accel);

axis_span = 2.5;

a = -axis_span;
b = axis_span;
ds = 1000;

x_ = linspace(-axis_span, axis_span, ds);
y_ = linspace(-axis_span, axis_span, ds);

[X,Y] = meshgrid(x_, y_);

for i = 1:ds
    for j = 1:ds
        x = [X(i, j), Y(i, j)];
        Z(i, j) = oracle(x);
    end
end

% Contour plot
hfig = my_figure();
contourf(X, Y, Z)
hold on;

curr = xs(2:end, :); 
prev = xs(1:end-1, :);

uv = curr - prev;
u = [uv(:, 1); 0];
v = [uv(:, 2); 0];

axis_lims = [a b a b];
hold on
plot(xs(:, 1), xs(:, 2), 'ko');
hold on
props = quiver(xs(:, 1), xs(:, 2), u, v, ...
               'color', 'green', 'AutoScale','off');
props.LineWidth = 2;
axis(axis_lims)
plot(xs(1, 1), xs(1, 2), 'kD', ...
                         'MarkerSize', 12, ...
                         'MarkerFaceColor','green');
hold on;
plot(xs(end, 1), xs(end, 2), '-s', ...
                             'MarkerSize', 12, ...
                             'MarkerFaceColor','red');

hold on;
plot(x_0, y_0, '-p', 'MarkerFaceColor','red',...
                     'MarkerSize',15);
hold off;

axis square
axis(axis_lims);

titletxt = sprintf(['$\\nu$ = ', num2str(nu), ', ', ...
                    '$\\sigma$ = ', num2str(sigma'), ', ', ...
                    '$\\lambda$ = ', num2str(lambda), ', ' ...
                    '$\\lambda_z$ = ', num2str(lambda_z)]);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_s = my_plot(iters, zbars, plot_config);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\Delta \hat{x}_1$', '$\Delta \hat{x}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_s = my_plot(iters, deltas, plot_config);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\hat{x}_1$', '$\hat{x}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_s = my_plot(iters, xs, plot_config);

% Save folder
path = [pwd '/../imgs/'];
fname = ['multiple_points', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu)];
saveas(hfig, [path, fname], 'epsc')

fname = ['zbars', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu), ...
         sprintf('lambz%.2f', lambda_z)];
saveas(hfigs_s, [path, fname], 'epsc')            
