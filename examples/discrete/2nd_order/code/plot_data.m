% Contour plot
hfig = my_figure();

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

contourf(X, Y, Z)
hold on;

for i = 1:length(xhat_tests)
    xs = xhat_tests{i};
    plot(xs(:, 1), xs(:, 2), 'r.');
    hold on;
end

curr = xhat_mean(2:end, :); 
prev = xhat_mean(1:end-1, :);

% X-Y plot
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

plot_config.titles = {titletxt, ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'Amplitude', 'Amplitude'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';
plot_config.markers = {'-', '--'};

[m, ~] = size(x_mean);
iters = (1:m)';
xhat_mean = xhat_mean(2:end, :);

hfigs_xmean = my_plot(iters, xhat_mean, plot_config);

plot_config.titles = {titletxt, ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [2, 1];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
hfigs_zmean = my_plot(iters, zbar_mean, plot_config);

plot_config.titles = {[titletxt, ' - $f(x, y) := ', latex(func),'$']};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$f(x, y)$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(fs_mean);
hfigs_fs = my_plot(iters, fs_mean, plot_config);

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

