close all
clear all
clc

% Method hyperparameter
nu = 1;
iterations = 50;

init_val = 1;
axis_span = 2;

sigma = 1;
zeta = 0;
lambda = 1;

% Save folder
path = [pwd '/images/'];

a = -1;
b = 1;
n = 1000;

oracle = @(x) x(1)^2 + x(2)^2;

% Batch version
x01 = uniform(a, b, [n, 2]);
x1 = expbary(oracle, x01, nu);

% Recursive version
x02 = [1, 1];
[ms, x2, xs] = drecexpbary(oracle, x02, nu, sigma, zeta, lambda, ...
                           iterations);

a = -axis_span;
b = axis_span;
n = 100;

x_ = linspace(-axis_span, axis_span, n);
y_ = linspace(-axis_span, axis_span, n);
[X,Y] = meshgrid(x_, y_);

for i = 1:length(x_)
    for j = 1:length(y_)
        Z(i, j) = oracle([X(i, j), Y(i, j)]);
    end
end

% Recursive version
hfig_recursive = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contour(X, Y, Z)
hold on;
plot(x2(1), x2(2), 'x');
hold on
plot(xs(:, 1), xs(:, 2), 'o');
hold on
plot(xs(:, 1), xs(:, 2), '-');
hold off

titletxt = sprintf(['Recursive barycenter method - $n$ = ', ...
                     num2str(iterations), ', ', ...
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
axis([-1, 1, -1, 1])

opt_msg = ['\leftarrow ' sprintf('(%.3f, %.3f)', x2(1), x2(2))];
text(x2(1), x2(2), opt_msg, 'FontSize', 15);

% Batched version
hfig_batch = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

plot(x01(:, 1), x01(:, 2), '*');
hold on
plot(x1(1), x1(2), 'x');
hold on

titletxt = sprintf(['Batched barycenter method - $n$ = ', ...
                     num2str(iterations), ', ', '$\\nu$ = ', ...
                     num2str(nu)], sigma, zeta, lambda);
htitle = title(titletxt); 
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

opt_msg = ['\leftarrow ' sprintf('(%.3f, %.3f)', x1(1), x1(2))];
text(x1(1), x1(2), opt_msg, 'FontSize', 15);

axis square
axis([-1, 1, -1, 1])

% Mass term
hfig_mval = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

plot(ms);

titletxt = sprintf(['Mass term $m_n$ - $n$ = ', ...
                     num2str(iterations), ', ', '$\\nu$ = ', ...
                     num2str(nu)], sigma, zeta, lambda);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('Iterações', 'interpreter', 'latex');
ylabel('$m_n$', 'interpreter', 'latex');

axis square

fname_batch = ['discrete_batch', ...
               '_nu', num2str(nu, '%.0e')];

saveas(hfig_batch, [path, fname_batch], 'epsc');
           
fname_recursive = ['discrete_recursive', ...
                   '_nu', num2str(nu, '%.0e'), ...
                   '_sigma', num2str(sigma, '%.0e'), ...
                   '_lambda', num2str(lambda, '%.0e')];

saveas(hfig_recursive, [path, fname_recursive], 'epsc');

fname_m = ['discrete_recursive_m', ...
           '_nu', num2str(nu, '%.0e'), ...
           '_sigma', num2str(sigma, '%.0e'), ...
           '_lambda', num2str(lambda, '%.0e')];

saveas(hfig_mval, [path, fname_m], 'epsc');



