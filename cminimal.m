close all
clear all
clc

n = 2;

% Method hyperparameter
nu = 1;
sigma = 1;
Sigma = sigma*ones(n, 1);
lambda = 0.98;

% Time integral
time = 0:0.01:100;

% Recursive version
oracle = @(x) sum(x.^2);
m0 = 1;
x0 = ones(n, 1);
[x, xs] = crecexpbary(oracle, m0, x0, nu, Sigma, lambda, time);

a = -1;
b = 1;
n = 100;

x_ = linspace(-1, 1);
y_ = linspace(-1, 1);
[X,Y] = meshgrid(x_, y_);
Z = X.^2 + Y.^2;

hfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contour(X, Y, Z)
hold on;
plot(x0(1), x0(2), '*');
hold on
plot(xs(:, 1), xs(:, 2), '-');
hold off
axis square

titletxt = sprintf(['$\\nu$ = ', num2str(nu), ', ', ...
                    '$\\sigma$ = [', num2str(Sigma'), ']$^{T}$, ', ...
                    '$\\lambda$ = ', num2str(lambda)]);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

axis square
opt_msg = ['\leftarrow ' sprintf('(%.3f, %.3f)', x(1), x(2))];
text(x(1), x(2), opt_msg, 'FontSize', 15);

% Save folder
path = [pwd '/images/'];
fname = ['continuous_nu', num2str(nu, '%.0e'), ...
         '_sigma', num2str(sigma, '%.0e'), ...
         '_lambda', num2str(lambda, '%.0e')];

saveas(hfig, [path, fname], 'epsc')
