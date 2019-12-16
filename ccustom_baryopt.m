close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 1;
sigma = 1;
lambda = 0;
zeta = 0;

% Time integral
time = 0:0.01:100;

% Recursive version
% oracle = @(x) (x(1) - 1)^2 + x(2)^2;
oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
              - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
              - 1/3*exp(-(x(1)+1).^2 - x(2).^2);
m0 = 1;
x0 = init_val*ones(n, 1);

curious_fun = @(t, xhat, m) curiosity_fun(t, xhat, m, oracle, ...
                                          nu, sigma, zeta);
[x, xs] = crecexpbary_custom(oracle, m0, x0, nu, lambda, ...
                             curious_fun, time);

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


hfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
contourf(X, Y, Z)
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
