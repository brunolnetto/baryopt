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
lambda_z = 0.5;
zeta = 2;

% Recursive version
syms x y
x_0 = 1;
y_0 = -1;
func = (x - x_0)^2 + (y - y_0)^2;
grad_f = gradient(func);
oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

FontSize = 25;

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

axis(axs{1}{1}, 'square');
axs{1}{1}.FontSize = FontSize;

delta_fs = [];
for i = 2:iterations
    delta_f_i = fs(i) - fs(i-1);
    delta_fs = [delta_fs; delta_f_i];
end

x0 = init_val*ones(n, 1);
m0 = 1;

iterations = 100;
n_iterations = 200;

is_accel = false;
options = struct('verbose', true);

xstar_tests = {};
xs_tests = {};
xhat_tests = {};
delta_tests = {};
zbars_tests = {}; 
F_tests = {};
ms_tests = {};
Fbar_tests = {};

accel_fun_txts = {'non_accel'};

wb = my_waitbar('Calculating minimum...');

for accel_fun_txt = accel_fun_txts
    accel_fun_txt = accel_fun_txt{1};
    
    for i = 1:n_iterations
        clear_inner_close_all(pwd);
        clear(accel_fun_txt);

        switch accel_fun_txt
            case 'pait_accel'
                accel_fun = @(m_1, x_1, delta_xhat_1) ...
                            pait_accel(m_1, x_1, delta_xhat_1, zeta);

            case 'integrated_accel'
                accel_fun = @(m_1, x_1, delta_xhat_1) ...
                              integrated_accel(m_1, x_1, delta_xhat_1, ...
                                               lambda_z, nu, oracle);
            case 'non_accel'
                accel_fun = @(m_1, x_1, delta_xhat_1) ...
                            non_accel(m_1, x_1, delta_xhat_1);

            otherwise
                warning('Accelerated non specified.');
                accel_fun = @(m_1, x_1, delta_xhat_1) ...
                            non_accel(m_1, x_1, delta_xhat_1);
        end

        options = struct('verbose', false, ...
                         'is_progress_visible', false);
        
        [x_star, xhats, xs, m, ...
         deltas, zbars, ...
         Fs_n, Fbars_n] = ...
         drecexpbary_custom(oracle, m0, x0, nu, sigma, ...
                            lambda, iterations, accel_fun, options);

        clear(func2str(accel_fun));

        xstar_tests{end+1} = x_star;
        xs_tests{end+1} = xs;
        xhat_tests{end+1} = xhats;
        delta_tests{end+1} = deltas;
        zbars_tests{end+1} = zbars;
        F_tests{end+1} = Fs_n;
        Fbar_tests{end+1} = Fbars_n;
        ms_tests{end+1} = m;

        wb = wb.update_waitbar(i, n_iterations);
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
    
    run('./plot_data.m');
end
