close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 1;
sigma = 1;
lambda = 0;

% Time integral
time = 0:0.01:40;

% Recursive version
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);
m0 = 1;
x0 = init_val*ones(n, 1);
zbar_0 = zeros(size(x0));

[x, xs, zs] = crecexpbary_custom(oracle, m0, x0, zbar_0, ...
                                 nu, lambda, sigma, time);

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


hfig_xy = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
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

hfig_z1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(time, zs(:, 1));

titletxt = 'Mean $\bar{z}_1$ of support density function';
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('t [s]', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');

hfig_z2 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(time, zs(:, 2));

titletxt = 'Mean $\bar{z}_2$ of support density function';
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('t [s]', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');

hfig_x1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(time, xs(:, 1));

titletxt = 'Barycenter $\hat{x}_1$ over time';
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('t [s]', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');

hfig_x2 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(time, xs(:, 2));

titletxt = 'Barycenter $\hat{x}_2$ over time';
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('t [s]', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');

hfig_m = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(time(1, 1:end-1), ms);

titletxt = 'Accumulated mass ';
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('t [s]', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');

% Save folder
path = [pwd '/images/'];
fname = ['continuous_nu', num2str(nu, '%.0e'), ...
         '_sigma', num2str(sigma, '%.0e'), ...
         '_lambda', num2str(lambda, '%.0e')];

saveas(hfig_xy, [path, fname], 'epsc')

path = [pwd '/images/'];
fname = 'curiosity1';

saveas(hfig_z1, [path, fname], 'epsc')

path = [pwd '/images/'];
fname = 'curiosity2';

saveas(hfig_z2, [path, fname], 'epsc')

path = [pwd '/images/'];
fname = 'barycenter1';

saveas(hfig_x1, [path, fname], 'epsc')

path = [pwd '/images/'];
fname = 'barycenter2';

saveas(hfig_x2, [path, fname], 'epsc')

path = [pwd '/images/'];
fname = 'mass';

saveas(hfig_m, [path, fname], 'epsc')
