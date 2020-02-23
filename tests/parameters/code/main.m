close all
clear all
clc

n = 2;

init_val = 0;

% Method hyperparameter
nu = 5;
sigma = 0.5;
lambda = 1;
lambda_z = 0.5;

% Time integral
iterations = 200;

% Recursive version
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

% syms x y
% func = (x - 1)^2 + (y + 1)^2;
% oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;
syms x y
func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
              - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
              - 1/3*exp(-(x+1)^2 - y^2);

oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
              - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
              - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = exp(-nu*oracle(x0));

n_iterations = 2000;

for nu = [1 3.3 10]
    for sigma = [0.1 1]
        for lambda = [0.9 1]
            x_test = [];

            wb = my_waitbar('Calculating minimum...');

            for i = 1:n_iterations
                [x, xs] = drecexpbary_custom(oracle, m0, x0, ...
                                             nu, sigma, lambda, ...
                                             lambda_z, iterations);
                x_test = [x_test, x];

                wb.update_waitbar(i, n_iterations);
            end

            axis_span = 2.5;

            a = -axis_span;
            b = axis_span;
            ds = 1000;

            x_ = linspace(-axis_span, axis_span, ds);
            y_ = linspace(-axis_span, axis_span, ds);

            [X,Y] = meshgrid(x_, y_);

            for i = 1:ds
                for j = 1:ds
                    x = [X(i, j), Y(i, j)];
                    Z(i, j) = oracle(x);
                end
            end

            hfig = my_figure();
            contourf(X, Y, Z)

            curr = xs(2:end, :); 
            prev = xs(1:end-1, :);

            uv = curr - prev;
            u = [uv(:, 1); 0];
            v = [uv(:, 2); 0];

            hold on
            plot(xs(:, 1), xs(:, 2), 'rx');
            hold on
            props = quiver(xs(:, 1), xs(:, 2), u, v, 'AutoScale','off');
            props.LineWidth = 2;
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

            ax = gca();
            tighten_plot(ax);

            % Save folder
            path = [pwd '/../imgs/'];
            fname = ['multiple_points', ...
                     sprintf('lamb%.2f', lambda), ...
                     sprintf('sigma%.2f', sigma), ...
                     sprintf('nu%.2f', nu)];
            fname = erase(fname,".");
            saveas(hfig, [path, fname, '.eps'], 'epsc')
        end
    end
end

