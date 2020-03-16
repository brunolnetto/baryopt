close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 0.5;
Sigma = sigma*eye(n);
lambda = 1;
lambda_z = 0.9;
zeta = 2;

% Recursive version
syms x y
x_0 = 1;
y_0 = -1;
func = (x - x_0)^2 + (y - y_0)^2;
grad_f = gradient(func);
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

% syms x y
% func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
%               - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
%               - 1/3*exp(-(x+1)^2 - y^2);
% 
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);
% 
% syms x y z w
% x_0 = 1;
% y_0 = 1;
% z_0 = 0;
% w_0 = 0;
% func = (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 + (w - w_0)^2;
% oracle =  @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2 + ...
%                (x(3) - z_0)^2 + (x(4) - w_0)^2;

x0 = init_val*ones(n, 1);
m0 = 1;

iterations = 100;
n_iterations = 100;

wb = my_waitbar('Calculating minimum...');
is_accel = false;

accel_fun_1 = @(m_1, xhat_1, delta_xhat_1) ...
                     non_accel(m_1, xhat_1, delta_xhat_1);
accel_fun_2 = @(m_1, xhat_1, delta_xhat_1) ...
            pait_accel(m_1, xhat_1, delta_xhat_1, zeta);
accel_fun_3 = @(m_1, xhat_1, delta_xhat_1) ...
                  integrated_accel(m_1, xhat_1, delta_xhat_1, ...
                                   lambda_z, nu, oracle);

accel_funs = {accel_fun_1, accel_fun_2, accel_fun_3};
accel_fun = accel_funs{1};

xstar_tests = {};
xs_tests = {};
xhat_tests = {};
delta_tests = {};
zbars_tests = {}; 
F_tests = {};
ms_tests = {};
Fbar_tests = {}; 
for i = 1:n_iterations            
    [x_star, xhats, xs, m, ...
     deltas, zbars, Fs_n, Fbars_n] = ...
     drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
                        lambda, iterations, accel_fun);
    
    clear(func2str(accel_fun));
                    
    xstar_tests{end+1} = x_star;
    xs_tests{end+1} = xs;
    xhat_tests{end+1} = xhats;
    delta_tests{end+1} = deltas;
    zbars_tests{end+1} = zbars;
    F_tests{end+1} = Fs_n;
    Fbar_tests{end+1} = Fbars_n;
    ms_tests{end+1} = m;
    
    wb.update_waitbar(i, n_iterations);
end

axis_span = 1.5;

a = -axis_span;
b = axis_span;
ds = 1000;

x_ = linspace(-axis_span, axis_span, ds);
y_ = linspace(-axis_span, axis_span, ds);

[X,Y] = meshgrid(x_, y_);

for i = 1:ds
    for j = 1:ds
        coords = [X(i, j), Y(i, j)];
        Z(i, j) = oracle(coords);
    end
end

% Contour plot
hfig = my_figure();
contourf(X, Y, Z)
hold on;

x_star = zeros(1, 2);
for i = 1:length(xhat_tests)
    x_aux = xstar_tests{i};
    x_star = x_star + x_aux;
end

x_star = x_star/n_iterations;

for i = 1:length(xhat_tests)
    xs = xhat_tests{i};
    plot(xs(:, 1), xs(:, 2), 'r.');
    hold on;
end

xhat_mean = zeros(size(xhat_tests{1}));
x_mean = zeros(size(xs_tests{1}));
F_mean = zeros(size(F_tests{1}));
Fbar_mean = zeros(size(Fbar_tests{1}));
for i = 1:length(xhat_tests)
    x_mean = x_mean + xs_tests{i};
    xhat_mean = xhat_mean + xhat_tests{i};
    F_mean = F_mean + F_tests{i};
    Fbar_mean = Fbar_mean + Fbar_tests{i};
end

xhat_mean = xhat_mean/n_iterations;
x_mean = x_mean/n_iterations;
F_mean = F_mean/n_iterations;
Fbar_mean = Fbar_mean/n_iterations;

curr = xhat_mean(2:end, :); 
prev = xhat_mean(1:end-1, :);

uv = curr - prev;
u = [uv(:, 1); 0];
v = [uv(:, 2); 0];

axis_lims = [a b a b];
hold on
plot(x_(:, 1), x_(:, 2), 'ko');
hold on
props = quiver(xhat_mean(:, 1), xhat_mean(:, 2), u, v, ...
               'color', 'green', 'AutoScale','off');
props.LineWidth = 2;
axis(axis_lims)
plot(x_star(1, 1), x_star(1, 2), ...
     'kD', 'MarkerSize', 12, ...
     'MarkerFaceColor','green');
hold on;
plot(x_star(1), x_star(2), 'kD',...
            'MarkerSize', 12, ...
            'MarkerFaceColor','green');
hold off;

axis square
axis(axis_lims);

titletxt = sprintf(['$\\nu$ = ', num2str(nu), ', ', ...
                    '$\\sigma$ = ', num2str(sigma'), ', ', ...
                    '$\\lambda$ = ', num2str(lambda), ', ' ...
                    '$\\lambda_z$ = ', num2str(lambda_z)]);
htitle = title(titletxt);
htitle.Interpreter = 'latex';
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');

ax = gca;
tighten_plot(ax);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'Amplitude', 'Amplitude'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'$\hat{x}_1$', '$x$'}, ...
                       {'$\hat{x}_2$', '$x$'}};
plot_config.pos_multiplots = [1, 2];
plot_config.markers = {{'-', '--'}, {'-', '--'}};

ys = {xhat_mean, x_mean};

iters = 1:length(xs(:, 1));
hfigs_xmean = my_plot(iters, xhat_mean, plot_config);

zbar_mean = zeros(size(zbars_tests{1}));
for i = 1:length(zbars_tests)
    zbar_mean = zbar_mean + zbars_tests{i};
end

zbar_mean = zbar_mean/length(zbars_tests);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_zmean = my_plot(iters, zbar_mean, plot_config);

fs = [];
for i = 1:iterations
    f_i = oracle(x_mean(i, :));
    fs = [fs; f_i];
end

plot_config.titles = {['$f(x, y) := ', latex(func),'$']};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$f(x, y)$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(fs);
hfigs_fs = my_plot(iters, fs, plot_config);


delta_mean = zeros(size(delta_tests{1}));
for i = 1:length(delta_tests)
    delta_mean = delta_mean + delta_tests{i};
end

delta_mean = delta_mean/n_iterations;
delta_mean = delta_mean(:, 1:2);

E_deltas = [];
[n_x, ~] = size(x_mean);
for i = 1:n_x
    grad_x_i = double(subs(grad_f, [x, y], x_mean(i, 1:2)));
    F_mean_i = F_mean(i);
    Fbar_mean_i = Fbar_mean(i);
    zbar_mean_i = zbar_mean(i, :);
    delta_xhat_i = F_mean_i*zbar_mean_i' - nu*Sigma*Fbar_mean_i*grad_x_i;
    
    E_deltas = [E_deltas; delta_xhat_i'];
end

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\Delta \hat{x}_1$', '$\Delta \hat{x}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.legends = {{'$\bar{\Delta \hat{x}_1}$', '$\rm{E}[\Delta \hat{x}_1]$'}, ...
                       {'$\bar{\Delta \hat{x}_2}$', '$\rm{E}[\Delta \hat{x}_2]$'}};
plot_config.pos_multiplots = [1, 2];
plot_config.markers = {{'-', '--'}, {'-', '--'}};

ys = {delta_mean, E_deltas};

iters = 1:length(delta_mean);
hfigs_delta = my_plot(iters, ys, plot_config);

% Save folder
path = [pwd '/../imgs/'];
posfix = [sprintf('lamb%d', 100*lambda), ...
         sprintf('sigma%d', 100*sigma), ...
         sprintf('nu%d', 100*nu), ...
         sprintf('lambz%d', 100*lambda_z)];
fname = ['multiple_points',  posfix];
saveas(hfig, [path, fname], 'epsc')

fname = ['x_mean', posfix];
saveas(hfigs_xmean, [path, fname], 'epsc')

fname = ['zbars', posfix];
saveas(hfigs_zmean, [path, fname], 'epsc')

fname = ['source', posfix];
saveas(hfigs_fs, [path, fname], 'epsc')

fname = ['deltas', posfix];
saveas(hfigs_delta, [path, fname], 'epsc')