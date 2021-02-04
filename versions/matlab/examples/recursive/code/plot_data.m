% Contour plot
hfig = my_figure();
contourf(X, Y, Z)
hold on;

n_sim = length(xhat_tests);
x_star = zeros(1, 2);
for i = 1:n_sim
    x_aux = xstar_tests{i};
    x_star = x_star + x_aux';
end

x_star = x_star/n_sim;

for i = 1:n_sim
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

xhat_mean = double(xhat_mean/n_sim);
x_mean = double(x_mean/n_sim);
F_mean = double(F_mean/n_sim);
Fbar_mean = double(Fbar_mean/n_iterations);

curr = xhat_mean(2:end, :); 
prev = xhat_mean(1:end-1, :);

uv = curr - prev;
u = [uv(:, 1); 0];
v = [uv(:, 2); 0];

FontSize = 20;

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

axs = gca;

axis(axs, 'square');
axs.FontSize = FontSize;
tighten_plot(axs);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'Amplitude', 'Amplitude'};
plot_config.grid_size = [1, 2];
plot_config.plot_type = 'stem';
plot_config.legends = {{'$\hat{x}_1$', '$x$'}, ...
                       {'$\hat{x}_2$', '$x$'}};
plot_config.pos_multiplots = [1, 2];
plot_config.markers = {{'-', '--'}, {'-', '--'}};

ys = {xhat_mean(1:end-1, :), x_mean};

iters = (1:(length(xs(:, 1)) - 1))';
[hfigs_xmean, axs] = my_plot(iters, ys, plot_config);

axis(axs{1}{1}, 'square');
axis(axs{1}{2}, 'square');
axs{1}{1}.FontSize = FontSize;
axs{1}{2}.FontSize = FontSize;

plot_config.titles = {''};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$F_n$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(F_mean);
[hfigs_Fmean, axs] = my_plot(iters, F_mean, plot_config);

axis(axs{1}{1}, 'square');
axs{1}{1}.FontSize = FontSize;

plot_config.titles = {''};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$\bar{F}_n$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(Fbar_mean);
[hfigs_Fbarmean, axs] = my_plot(iters, Fbar_mean, plot_config);

axis(axs{1}{1}, 'square');
axs{1}{1}.FontSize = FontSize;

zbar_mean = zeros(size(zbars_tests{1}));
for i = 1:length(zbars_tests)
    zbar_mean = zbar_mean + zbars_tests{i};
end

zbar_mean = zbar_mean/length(zbars_tests);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
plot_config.grid_size = [1, 2];
plot_config.plot_type = 'stem';

iters = 1:length(zbars(:, 1));
[hfigs_zmean, axs] = my_plot(iters, zbar_mean, plot_config);

axis(axs{1}{1}, 'square');
axis(axs{1}{2}, 'square');
axs{1}{1}.FontSize = FontSize;
axs{1}{2}.FontSize = FontSize;

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
[hfig_fvalue, axs] = my_plot(iters, fs, plot_config);

axis(axs{1}{1}, 'square');
axs{1}{1}.FontSize = FontSize;

delta_fs = [];
for i = 2:iterations
    delta_f_i = fs(i) - fs(i-1);
    delta_fs = [delta_fs; delta_f_i];
end

plot_config.titles = {['']};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$\Delta \, f(x, y)$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(fs);
[hfigs_fs, axs] = my_plot(iters(2:end), delta_fs, plot_config);

xhat_0 = [0; 0];
xhat_1 = xhat_0;
m_1 = 0;

xhats = [];

for i = 1:iterations
    grad_f_val = my_subs(grad_f, [x; y], xhat_1);
    f_val = my_subs(func, [x; y], xhat_1);
    
    e_i = exp(-nu*f_val);
    m = m_1 + e_i;
    
    F_n = e_i/m;
    Fbar_n = (m_1/m)*F_n;
    
    xhat = xhat_1 - nu*sigma*Fbar_n*grad_f_val;
    
    xhat_1 = xhat;
    m_1 = m;
    
    xhats = [xhats; xhat'];
end

fvals_determ = [];
for i = 1:iterations
    f_val = double(subs(func, [x; y], xhats(i, :)'));
    fvals_determ = [fvals_determ; f_val];
end

plot_config.titles = {['$f(x, y) := ', latex(func),'$']};
plot_config.xlabels = {'Iterations'};
plot_config.ylabels = {'$f(x, y)$'};
plot_config.grid_size = [1, 1];
plot_config.plot_type = 'stem';

iters = 1:length(fs);
[hfig_fvalue_det, axs] = my_plot(iters, fvals_determ, plot_config);


axis(axs{1}{1}, 'square');
axs{1}{1}.FontSize = FontSize;

delta_mean = zeros(size(delta_tests{1}));
for i = 1:length(delta_tests)
    delta_mean = delta_mean + delta_tests{i};
end

delta_mean = delta_mean/n_iterations;
delta_mean = delta_mean(:, 1:2);

plot_config.titles = {'', ''};
plot_config.xlabels = {'', 'Iterations'};
plot_config.ylabels = {'$\Delta \hat{x}_1$', '$\Delta \hat{x}_2$'};
plot_config.grid_size = [1, 2];
plot_config.plot_type = 'stem';

iters = 1:length(delta_mean(:, 1));
[hfigs_deltamean, axs] = my_plot(iters, delta_mean, plot_config);

axis(axs{1}{1}, 'square');
axis(axs{1}{2}, 'square');
axs{1}{1}.FontSize = FontSize;
axs{1}{2}.FontSize = FontSize;

% Save folder
path = [pwd '/../imgs/'];

if(~contains(path, accel_fun_txt))
    mkdir([path, accel_fun_txt])
end

posfix = [sprintf('lamb%d', 100*lambda), ...
          sprintf('sigma%d', 100*sigma), ...
         sprintf('nu%d', 100*nu), ...
         sprintf('lambz%d', 100*lambda_z)];
fname = ['multiple_points',  posfix];
saveas(hfig, [path, fname], 'epsc')

fname = [accel_fun_txt, '/', 'f_value', accel_fun_txt, '_', posfix];
saveas(hfig_fvalue_det, [path, fname], 'epsc')

fname = [accel_fun_txt, '/', 'f_value', accel_fun_txt, '_', posfix];
saveas(hfig_fvalue, [path, fname], 'epsc')

fname = [accel_fun_txt, '/', 'x_mean', accel_fun_txt, '_', posfix];
saveas(hfigs_xmean, [path, fname], 'epsc')

fname = [accel_fun_txt, '/', 'zbars', accel_fun_txt, '_', posfix];
saveas(hfigs_zmean, [path, fname], 'epsc')

fname = [accel_fun_txt, '/','delta', accel_fun_txt, '_', posfix];
saveas(hfigs_deltamean, [path, fname], 'epsc')

fname = [accel_fun_txt, '/', 'source', accel_fun_txt, '_', posfix];
saveas(hfigs_fs, [path, fname], 'epsc')
