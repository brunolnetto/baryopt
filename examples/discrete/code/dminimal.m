close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 1;
lambda = 1;
lambda_z = 0;

% Recursive version
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

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
n_iterations = 500;

wb = my_waitbar('Calculating minimum...');
is_accel = false;

x_tests = {};
for i = 1:n_iterations
    [x, xs, m, ...
     deltas, zbars] = drecexpbary_custom(oracle, m0, x0, ...
                                         nu, sigma, lambda, ...
                                         lambda_z, iterations, ...
                                         is_accel);
    x_tests{end+1} = xs;
    
    wb.update_waitbar(i, n_iterations);
end

axis_span = 1.5;

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

x_mean = zeros(1, 2);
for i = 1:length(x_tests)
    x = x_tests{i};
    x_mean = x_mean + x(end, :);
end

x_mean = x_mean/n_iterations;

for i = 1:length(x_tests)
    xs = x_tests{i};
    plot(xs(:, 1), xs(:, 2), 'r.');
    hold on;
end

x_ = zeros(size(x_tests{1}));
for i = 1:length(x_tests)
    x_ = x_ + x_tests{i};
end

x_ = x_/n_iterations;

curr = x_(2:end, :); 
prev = x_(1:end-1, :);

uv = curr - prev;
u = [uv(:, 1); 0];
v = [uv(:, 2); 0];

axis_lims = [a b a b];
hold on
plot(x_(:, 1), x_(:, 2), 'ko');
hold on
props = quiver(x_(:, 1), x_(:, 2), u, v, ...
               'color', 'green', 'AutoScale','off');
props.LineWidth = 2;
axis(axis_lims)
plot(x_(1, 1), x_(1, 2), 'kD', 'MarkerSize', 12, ...
                         'MarkerFaceColor','green');
hold on;
plot(x_mean(1), x_mean(2), 'kD',...
            'MarkerSize', 12, ...
            'MarkerFaceColor','green');
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

ax = gca;
tighten_plot(ax);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_s = my_plot(iters, zbars, plot_config);

% Save folder
path = [pwd '/imgs/'];
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
