close all
clear all
clc

n = 2;

init_val = 1;

% Method hyperparameter
nu = 0.5;
sigma = 1;
lambda = 1;

% Time integral
iterations = 2000;

% Recursive version
% oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;
oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
              - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
              - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = exp(-nu*oracle(x0));

[x, xs] = drecexpbary_custom(oracle, m0, x0, nu, sigma, lambda, iterations);

axis_span = 2;

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
contour(X, Y, Z)
hold on;
plot(x0(1), x0(2), '*');
hold on
plot(xs(:, 1), xs(:, 2), 'o-');
hold off
axis square
axis([a b a b])

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
fname = 'continuous_nu';

saveas(hfig, [path, fname], 'epsc')
