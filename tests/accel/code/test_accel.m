close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 0.5;
lambda = 1;

% Recursive version
x_0 = 1;
y_0 = -1;
oracle = @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;

syms x y
func = (x - 1)^2 + (y + 1)^2;
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;
% syms x y
% func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
%               - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
%               - 1/3*exp(-(x+1)^2 - y^2);
% 
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = exp(-nu*oracle(x0));

iterations = 20;
n_iterations = 1000;

x_test = {};

wb = my_waitbar('Calculating minimum...');
is_accels = [true];

for is_accel = is_accels
    if(~is_accel)
        lambda_z = 0;
    else
        lambda_z = 0.9;
    end
    
    zbars_total = zeros(iterations, 2);
    for i = 1:n_iterations
        [x, xs] = drecexpbary_custom(oracle, m0, x0, ...
                                     nu, sigma, lambda, lambda_z, ...
                                     iterations, is_accel);
        x_test{end+1} = xs;
        
        zbars_total = zbars_total + zbars;
        
        wb.update_waitbar(i, n_iterations);
    end
    
    zbars = zbars_total/iterations;
    
    axis_span = 1.5;
    a = -axis_span;
    b = axis_span;
    ds = 1000;

    axis_lims = [a+x_0 b+x_0 a+y_0 b+y_0];
    x_ = linspace(axis_lims(1), axis_lims(2), ds);
    y_ = linspace(axis_lims(3), axis_lims(4), ds);

    [X,Y] = meshgrid(x_, y_);

    for i = 1:ds
        for j = 1:ds
            x = [X(i, j), Y(i, j)];
            Z(i, j) = oracle(x);
        end
    end

    % Contour plot
    hfig = my_figure();
    contourf(X, Y, Z)
    hold on;
    
    x_mean = zeros(1, 2);
    for i = 1:length(x_test)
        x = x_test{i};
        x_mean = x_mean + x(end, :);
    end

    x_mean = x_mean/n_iterations;

    for i = 1:length(x_test)
        xs = x_test{i};
        plot(xs(:, 1), xs(:, 2), 'r.');
        hold on;
    end
    
    x_ = zeros(size(x_test{1}));
    for i = 1:length(x_test)
        x_ = x_ + x_test{i};
    end

    x_ = x_/n_iterations;

    curr = x_(2:end, :); 
    prev = x_(1:end-1, :);

    uv = curr - prev;
    u = [uv(:, 1); 0];
    v = [uv(:, 2); 0];

    hold on
    plot(x_(:, 1), x_(:, 2), 'ko');
    hold on
    props = quiver(x_(:, 1), x_(:, 2), u, v, ...
                   'color', 'green', 'AutoScale','off');
    props.LineWidth = 2;
    axis(axis_lims)
    plot(x_(1, 1), x_(1, 2), 'kD', 'MarkerSize', 12, ...
                             'MarkerFaceColor','green');
    hold on;
    plot(x_mean(1), x_mean(2), 'kD',...
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
    
    % Acceleration plot
    plot_config.titles = {'', ''};
    plot_config.xlabels = {'', 'Iterations'};
    plot_config.ylabels = {'$\bar{z}_1$', '$\bar{z}_2$'};
    plot_config.grid_size = [2, 1];
    plot_config.plot_type = 'stem';

    iters = 1:length(zbars(:, 1));
    hfigs_s = my_plot(iters, zbars, plot_config);

    % Save folder
    path = [pwd '/../imgs/'];
    fname = ['multiple_points_', ...
         sprintf('lamb%.2f', lambda), ...
         sprintf('sigma%.2f', sigma), ...
         sprintf('nu%.2f', nu)];
    saveas(hfig, [path, fname], 'epsc')

    fname = ['zbars_', ...
             sprintf('lamb%.2f', lambda), ...
             sprintf('sigma%.2f', sigma), ...
             sprintf('nu%.2f', nu), ...
             sprintf('lambz%.2f', lambda_z)];
    saveas(hfigs_s, [path, fname], 'epsc')
end