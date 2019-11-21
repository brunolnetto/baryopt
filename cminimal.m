close all
clear all
clc

% Method hyperparameter
nu = 1;
Sigma = [0.2; 0.2];
lambda = 1;

% Time integral
time = 0:0.1:1;

% Save folder
path = [pwd '/images/'];
fname = 'continuous';

a = -1;
b = 1;
n = 100;

% Recursive version
oracle = @(x) sum(x.^2);
x0 = [1; 1];
[x, xs] = crecexpbary(oracle, x0, nu, Sigma, lambda, time);

x = linspace(-1, 1);
y = linspace(-1, 1);
[X,Y] = meshgrid(x,y);
Z = X.^2 + Y.^2;

hfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contour(X, Y, Z)
hold on;
plot(x01(:, 1), x01(:, 2), '*');
hold on
plot(x1(1), x1(2), 'x');
hold on
plot(xs(:, 1), xs(:, 2));
hold on
plot(xs(:, 1), xs(:, 2), 'o');
hold off

titletxt = sprintf(['$n$ = ', num2str(iterations), ', ', ...
                    '$\\nu$ = ', num2str(nu), ', ', ...
                    '$\\sigma$ = ', num2str(sigma), ', ', ...
                    '$\\xi$ = ', num2str(zeta), ', ', ...
                    '$\\lambda$ = ', num2str(lambda)], ...
                    sigma, zeta, lambda);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

axis square

saveas(hfig, [path, fname], 'epsc')
