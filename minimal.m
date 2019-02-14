close all
clear all
clc

oracle = @(x) sum(x.^2);

nu = 2;
iterations = 10000;

sigma = 0.5;
zeta = 0;
lambda = 0.95;
gamma = 0.5;

path = [pwd '/images/'];
fname = 'shapefactor';

a = -10;
b = 10;
n = 100;

% Batch version
x01 = uniform(a, b, [n, 2]);
x1 = expbary(oracle, x01, nu);

% Recursive version
oracle = @(x) sum(x.^2);
x02 = [10, 10];
[x2, xs] = recexpbary(oracle, x02, nu, ...
                      sigma, zeta, lambda, gamma, ...
                      iterations, 'shape');

x = linspace(-10, 10);
y = linspace(-10, 10);
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

titletxt = sprintf(['$\\sigma$ = ', num2str(sigma), ', ', ...
                    '$\\xi$ =', num2str(zeta), ', ', ...
                    '$\\lambda$ =', num2str(lambda), ', ', ...
                    '$\\gamma$ =', num2str(gamma)], sigma, zeta, lambda, gamma);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x');
ylabel('y');

axis square

saveas(hfig, [path, fname], 'epsc')
