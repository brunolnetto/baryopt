oracle = @(x) sum(x.^2);

nu = 1;
error = 1e-3;
iterations = 10000;
sigma = 0.1;
zeta = 0.1;
lambda = 0.5;

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
                      sigma, zeta, lambda, iterations);

x = linspace(-10, 10);
y = linspace(-10, 10);
[X,Y] = meshgrid(x,y);
Z = X.^2 + Y.^2;

figure
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

axis square