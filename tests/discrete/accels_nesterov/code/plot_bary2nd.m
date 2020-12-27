% -------- Plots --------
title_mask = '$\\nu = %.2f$, $\\sigma = %.2f$, $\\zeta = %.2f$, $\\lambda_z = %.2f$';
title = sprintf(title_mask, nu, sigma, zeta, lambda_z);

legend_atom = {'$\bar{z} = 0$', ...
               '$\bar{z} = \zeta \, \Delta \hat{x}$', ...
               '$\bar{z} = barycenter([\Delta \hat{x}_1, \cdots, \Delta \hat{x}_{n-1}])$'};

markers_atom = {'-', '--', '.-'};

% zbar plot
plot_config_z.titles = {'', ''};
plot_config_z.xlabels = {'', 'Iterations'};
plot_config_z.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config_z.grid_size = [2, 1];
plot_config_z.plot_type = 'stem';
plot_config_z.legends = {legend_atom, legend_atom};
plot_config_z.pos_multiplots = [1, 1, 2, 2];
plot_config_z.markers = {markers_atom, markers_atom};

zbars_means23 = double([zbars_means{2}, zbars_means{3}]);
zbars_ = {double(zbars_means{1}), ...
          [zbars_means23(:, 1), zbars_means23(:, 3), ...
           zbars_means23(:, 2), zbars_means23(:, 4)]};

[m, n] = size(zbars(:, 1));
iters = (1:length(zbars(:, 1)))';
hfigs_zmean = my_plot(iters, zbars_, plot_config_z);

% function plot
fs_ = {double(fs_means{1}), ...
       double([fs_means{2}, fs_means{3}])};

plot_config_f.titles = {[title, ' $f(x, y) := ', latex(func),'$']};
plot_config_f.xlabels = {'Iterations'};
plot_config_f.ylabels = {'$f(x, y)$'};
plot_config_f.grid_size = [1, 1];
plot_config_f.legends = legend_atom;
plot_config_f.pos_multiplots = [1, 1];
plot_config_f.markers = markers_atom;

[m, n] = size(fs_{1});
iters = 1:m;
hfigs_fs = my_plot(iters, fs_, plot_config_f);

% Delta xhat plot
deltas_means23 = double([deltas_means{2}, deltas_means{3}]);
delta_xhat_ = {double(deltas_means{1}), ...
               double([deltas_means23(:, 1), deltas_means23(:, 3), ...
                       deltas_means23(:, 2), deltas_means23(:, 4)])};

plot_config_delta.titles = {title, ''};
plot_config_delta.xlabels = {'', 'Iterations'};
plot_config_delta.ylabels = {'$\Delta \hat{x}_1$', '$\Delta \hat{x}_2$'};
plot_config_delta.grid_size = [2, 1];
plot_config_delta.plot_type = 'stem';
plot_config_delta.legends = {legend_atom, legend_atom};
plot_config_delta.pos_multiplots = [1, 1, 2, 2];
plot_config_delta.markers = {markers_atom, markers_atom};

iters = (1:length(delta_xhat_{1}))';
hfigs_delta_mean = my_plot(iters, delta_xhat_, plot_config_delta);

% xhat plot
xhat_means23 = double([xhat_means{2}, xhat_means{3}]);
xhat_ = {double(xhat_means{1}), ...
          double([xhat_means23(:, 1), xhat_means23(:, 3), ...
                  xhat_means23(:, 2), xhat_means23(:, 4)])};

plot_config_xhat.titles = {title, ''};
plot_config_xhat.xlabels = {'', 'Iterations'};
plot_config_xhat.ylabels = {'$\hat{x}_1$', '$\hat{x}_2$'};
plot_config_xhat.grid_size = [2, 1];
plot_config_xhat.plot_type = 'stem';
plot_config_xhat.legends = {legend_atom, legend_atom};
plot_config_xhat.pos_multiplots = [1, 1, 2, 2];
plot_config_xhat.markers = {markers_atom, markers_atom};

iters = (1:length(xhat_{1}))';
hfigs_xhatmean = my_plot(iters, xhat_, plot_config_xhat);

% x plot
x_means23 = double([x_means{2}, x_means{3}]);
x_ = {double(x_means{1}), double([x_means23(:, 1), x_means23(:, 3), ...
                          x_means23(:, 2), x_means23(:, 4)])};

plot_config_xhat.titles = {title, ''};
plot_config_xhat.xlabels = {'', 'Iterations'};
plot_config_xhat.ylabels = {'$x_1$', '$x_2$'};
plot_config_xhat.grid_size = [2, 1];
plot_config_xhat.plot_type = 'stem';
plot_config_xhat.legends = {legend_atom, legend_atom};
plot_config_xhat.pos_multiplots = [1, 1, 2, 2];
plot_config_xhat.markers = {markers_atom, markers_atom};

iters = (1:length(x_{1}))';
hfigs_xmean = my_plot(iters, x_, plot_config_xhat);

% y plot
y_means23 = double([y_means{2}, y_means{3}]);
y_ = {double(y_means{1}), double([y_means23(:, 1), y_means23(:, 3), ...
                                  y_means23(:, 2), y_means23(:, 4)])};

plot_config_xhat.titles = {title, ''};
plot_config_xhat.xlabels = {'', 'Iterations'};
plot_config_xhat.ylabels = {'$y_1$', '$y_2$'};
plot_config_xhat.grid_size = [2, 1];
plot_config_xhat.plot_type = 'stem';
plot_config_xhat.legends = {legend_atom, legend_atom};
plot_config_xhat.pos_multiplots = [1, 1, 2, 2];
plot_config_xhat.markers = {markers_atom, markers_atom};

iters = (1:length(y_{1}))';
hfigs_ymean = my_plot(iters, y_, plot_config_xhat);

% m plot
m_ = {double(ms_means{1}), double([ms_means{2}, ms_means{3}])};

plot_config_m.titles = {''};
plot_config_m.xlabels = {'Iterations'};
plot_config_m.ylabels = {'$x_1$'};
plot_config_m.grid_size = [1, 1];
plot_config_m.plot_type = 'stem';
plot_config_m.legends = {legend_atom};
plot_config_m.pos_multiplots = [1, 1];
plot_config_m.markers = {markers_atom};

iters = (1:length(m_{1}))';
hfigs_mmean = my_plot(iters, m_, plot_config_m);

% Save folder
path = [pwd '/../imgs/'];
posfix = [sprintf('lamb%d', 100*lambda), ...
          sprintf('sigma%d', 100*sigma), ...
          sprintf('nu%d', 100*nu), ...
          sprintf('lambz%d', 100*lambda_z)];

fname = ['xmean_', posfix];
saveas(hfigs_xmean, [path, fname], 'epsc')

fname = ['ymean_', posfix];
saveas(hfigs_ymean, [path, fname], 'epsc')

fname = ['xhatmean_', posfix];
saveas(hfigs_xhatmean, [path, fname], 'epsc')

fname = ['zbars_', posfix];
saveas(hfigs_zmean, [path, fname], 'epsc')

fname = ['source_', posfix];
saveas(hfigs_fs, [path, fname], 'epsc')

fname = ['ms_', posfix];
saveas(hfigs_mmean, [path, fname], 'epsc')