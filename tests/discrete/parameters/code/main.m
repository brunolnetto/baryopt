close all
clear all
clc

n = 2;

init_val = 0;

% Time integral
iterations = 50;

% Recursive version
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

syms x y
func = (x - 1)^2 + (y + 1)^2;
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

% syms x y
% func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
%        - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
%        - 1/3*exp(-(x+1)^2 - y^2);
% 
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = 0;

n_iterations = 75;

accel_fun = @(m_1, xhat_1, ...
              delta_xhat_1) ...
              non_accel(m_1, xhat_1, delta_xhat_1);

nus = [1, 10, 20];
sigmas = [0.1, 1];
lambdas = [1, 0.5];

j = 1;
len_params = length(nus)*length(sigmas)*length(lambdas);

great_wb = my_waitbar('Parameter great waitbar');

for nu = nus
    for sigma = sigmas
        for lambda = lambdas
            xhat_tests = {};

            wb = my_waitbar('Calculating minimum...');

            for i = 1:n_iterations
                [x_star, xhats, xs, ms, ...
                 deltas, zbars, Fs_n, Fbars_n] = ...
                    drecexpbary_custom(oracle, m0, x0, ...
                                       nu, sigma, lambda, ...
                                       iterations, accel_fun);
                
                xhat_tests{i} = xhats;

                wb.update_waitbar(i, n_iterations);
            end
            
            wb.close_window();
            
            % Average value
            xhat_mean = zeros(size(xhat_tests{1}));
            for i = 1:length(xhat_tests)
                xhat_mean = xhat_mean + xhat_tests{i};
            end
            xhat_mean = xhat_mean/n_iterations;
            
            axis_span = 3.5;

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

            curr = xhat_mean(2:end, :); 
            prev = xhat_mean(1:end-1, :);

            uv = curr - prev;
            u = [uv(:, 1); 0];
            v = [uv(:, 2); 0];
            
            for i = 1:length(xhat_tests)
                hold on
                xhat_test = xhat_tests{i};
                plot(xhat_test(:, 1), xhat_test(:, 2), 'rx');
            end
            
            hold on
            plot(xhat_mean(:, 1), xhat_mean(:, 2), 'rx');
            hold on
            props = quiver(xhat_mean(:, 1), xhat_mean(:, 2), ...
                           u, v, 'AutoScale','off');
            props.LineWidth = 2;
            hold on
            plot(xhat_mean(end, 1), xhat_mean(end, 2), ...
                 '-p', 'MarkerFaceColor', 'red', 'MarkerSize',15);
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
            ax.Units = hfig.Units;
            ax.FontSize = 25;
            
            % Save folder
            path = [pwd '/../imgs/'];
            fname = ['multiple_points', ...
                     sprintf('lamb%.2f', lambda), ...
                     sprintf('sigma%.2f', sigma), ...
                     sprintf('nu%.2f', nu)];
            fname = erase(fname,".");
            saveas(hfig, [path, fname, '.eps'], 'epsc')
            saveas(hfig, [path, fname, '.fig'])
            
            great_wb.update_waitbar(j, len_params);
            j = j + 1;
        end
    end
end

great_wb.close_window();
