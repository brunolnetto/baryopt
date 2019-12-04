close all
clear all
clc

n = 2;

init_val = 1;
axis_span = 2;

% Method hyperparameter
nu = 1;
sigma = 1;
omega = 1;
lambda = 0;

% Time integral
time = 0:0.001:10;

% Recursive version
% oracle = @(x) sum(x.^2);
oracle = @(x) (x(2))^2;
m0 = 1;
x0 = init_val*ones(n, 1);

curious_fun = @(t, xhat) curiosity_fun(t, xhat, sigma);
[x, xs] = crecexpbary_custom(oracle, m0, x0, nu, lambda, ...
                             curious_fun, time);

a = -axis_span;
b = axis_span;
n = 100;

x_ = linspace(-axis_span, axis_span, n);
y_ = linspace(-axis_span, axis_span, n);
[X,Y] = meshgrid(x_, y_);

for i = 1:length(x_)
    for j = 1:length(y_)
        Z(i, j) = oracle([x_(i), y_(j)]);
    end
end


hfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contour(X, Y, Z)
hold on;
plot(x0(1), x0(2), '*');
hold on
plot(xs(:, 1), xs(:, 2), '-');
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
opt_msg = ['  \leftarrow ' sprintf('(%.2f, %.2f)', x(1), x(2))];
text(x(1), x(2), opt_msg, 'FontSize', 15);

% Save folder
path = [pwd '/images/'];
fname = ['continuous_nu', num2str(nu, '%.0e'), ...
         '_sigma', num2str(sigma, '%.0e'), ...
         '_lambda', num2str(lambda, '%.0e')];

saveas(hfig, [path, fname], 'epsc')
