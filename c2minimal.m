close all
clear all
clc

n = 2;

init_val = 1;

% Method hyperparameter
nu = 0.01;
sigma = 1;
lambda = 1;

% Time integral
iterations = 1000;

% Recursive version
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;
m0 = 1;
x0 = init_val*ones(n, 1);

[x, xs] = d2recexpbary_custom(oracle, m0, x0, ...
                              nu, sigma, lambda, iterations);

axis_span = 5;

a = -axis_span;
b = axis_span;
ds = 100;

x_ = linspace(-axis_span, axis_span, ds);
y_ = linspace(-axis_span, axis_span, ds);

[X,Y] = meshgrid(x_, y_);

for i = 1:ds
    for j = 1:ds
        x = [X(i, j), Y(i, j)];
        Z(i, j) = oracle(x);
    end
end

hfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contourf(X, Y, Z)
hold on;
plot(x0(1), x0(2), '*');
hold on
plot(xs(:, 1), xs(:, 2), 'o-');
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
opt_msg = ['\leftarrow ' sprintf('(%.3f, %.3f)', xs(end, 1), xs(end, 2))];
text(xs(end, 1), xs(end, 2), opt_msg, 'FontSize', 15);

% Save folder
path = [pwd '/images/'];
fname = ['continuous_nu', num2str(nu, '%.0e'), ...
         '_sigma', num2str(sigma, '%.0e'), ...
         '_lambda', num2str(lambda, '%.0e')];

saveas(hfig, [path, fname], 'epsc')
