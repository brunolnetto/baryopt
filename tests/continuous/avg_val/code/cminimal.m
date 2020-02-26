close all
clear all
clc

n = 2;

init_val = 0;

% Recursive version
x_0 = 1;
y_0 = -1;
oracle = @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

n_iterations = 50;

% Method hyperparameter
nu = 10;
sigma = 1;
lambda = 0;
lambda_z = 0;

% Time integral
time = 0:0.2:200;
time = time';

m0 = 1;
x0 = init_val*ones(n, 1);
zbar_0 = zeros(size(x0));

x_tests = {};
z_tests = {};
m_tests = {};

for i = 1:n_iterations
    [x, xs, zs, ms] = ...
            crecexpbary_custom(oracle, m0, x0, zbar_0, nu, ...
                               lambda, lambda_z, sigma, time);
    x_tests{end+1} = xs;
    z_tests{end+1} = zs;
    m_tests{end+1} = ms;
end

x_mean = zeros(size(x_tests{1}));
z_mean = zeros(size(z_tests{1}));
m_mean = zeros(size(m_tests{1}));
for i = 1:n_iterations
    x_mean = x_mean +  x_tests{i};
    z_mean = z_mean +  z_tests{i};
    m_mean = m_mean +  m_tests{i};
end

x_mean = x_mean/n_iterations;
z_mean = z_mean/n_iterations;
m_mean = m_mean/n_iterations;

m_mean = m_mean';

axis_span = 2;
a = -axis_span;
b = axis_span;
n = 100;

x_ = linspace(-axis_span, axis_span, n);
y_ = linspace(-axis_span, axis_span, n);
[X,Y] = meshgrid(x_, y_);

[m, n] = size(X);

for i = 1:m
    for j = 1:n
        Z(i, j) = oracle([X(i, j), Y(i, j)]);
    end
end

hfig_xy = my_figure();
contourf(X, Y, Z)
hold on;
plot(x0(1), x0(2), ...
     '-s', 'MarkerFaceColor','green',...
     'MarkerSize',15);
hold on
plot(x_mean(:, 1), x_mean(:, 2), '-');
hold on
plot(x_0, y_0, '-p', 'MarkerFaceColor','red', 'MarkerSize',15);
hold off
axis square

titletxt = sprintf(['$\\nu$ = ', num2str(nu), ', ', ...
                    '$\\sigma$ = ', num2str(sigma'), ', ', ...
                    '$\\lambda$ = ', num2str(lambda)]);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

axis square
opt_msg = ['  \leftarrow ' sprintf('(%.2f, %.2f)', ...
           x_mean(end, 1), x_mean(end, 2)), ];
text(x_mean(end, 1), x_mean(end, 2), opt_msg, 'FontSize', 15);

% x-mean coordinates
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 't'};
plot_config.ylabels = {'$\bar{\hat{x}}_1$', '$\bar{\hat{x}}_2$'};
plot_config.grid_size = [2, 1];

hfig_x = my_plot(time, x_mean, plot_config);

% z-mean coordinates
plot_config.titles = {'', ''};
plot_config.xlabels = {'', 't'};
plot_config.ylabels = {'$\bar{\bar{z}}_1$', '$\bar{\bar{z}}_2$'};
plot_config.grid_size = [2, 1];

hfig_z = my_plot(time, z_mean, plot_config);

% m-mean coordinates
plot_config.titles = {''};
plot_config.xlabels = {'t'};
plot_config.ylabels = {'$\bar{m}$'};
plot_config.grid_size = [1, 1];

hfig_m = my_plot(time, m_mean, plot_config);
axis square;

% Save folder
path = [pwd '/imgs/'];
fname = 'xy';
saveas(hfig_xy, [path, fname], 'epsc')

fname = 'pos';
saveas(hfig_x, [path, fname], 'epsc')

fname = 'accel';
saveas(hfig_z, [path, fname], 'epsc')

fname = 'mass';
saveas(hfig_m, [path, fname], 'epsc')
