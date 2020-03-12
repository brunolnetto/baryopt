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
% syms x y
% x_0 = 1;
% y_0 = -1;
% func = (x - x_0)^2 + (y - y_0)^2;
% grad_f = gradient(func);
% oracle = @(x) (x(1) - 1)^2 + (x(2) + 1)^2;

syms x y
x_0 = 1;
y_0 = -1;
func = (x - x_0)^2 + (y - y_0)^2;
oracle =  @(x) (x(1) - x_0)^2 + (x(2) - y_0)^2;

% syms x y
% func = 3*(1-x).^2.*exp(-(x^2) - (y+1)^2) ... 
%               - 10*(x/5 - x.^3 - y^5).*exp(-x^2-y^2) ... 
%               - 1/3*exp(-(x+1)^2 - y^2);
% 
% oracle = @(x) 3*(1-x(1)).^2.*exp(-(x(1)^2) - (x(2)+1).^2) ... 
%               - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
%               - 1/3*exp(-(x(1)+1).^2 - x(2).^2);

x0 = init_val*ones(n, 1);
m0 = 1;

iterations = 200;
n_iterations = 1000;

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
accel_fun = accel_funs{3};

xstar_tests = {};
xs_tests = {};
ys_tests = {};
xhat_tests = {};
delta_tests = {};
zbars_tests = {}; 
ms_tests = {};
for i = 1:n_iterations
    [x_star, xhats, xs, ys, ms, deltas, zbars] = ...
     drecexpbary_2nd_custom(oracle, m0, x0, nu, sigma, ...
                            lambda, iterations, accel_fun);
    
    clear(func2str(@integrated_accel));
    
    xstar_tests{end+1} = x_star;
    xs_tests{end+1} = xs;
    ys_tests{end+1} = ys;
    xhat_tests{end+1} = xhats;
    delta_tests{end+1} = deltas;
    zbars_tests{end+1} = zbars;    
    ms_tests{end+1} = ms;
    
    wb.update_waitbar(i, n_iterations);
end

% run('./average_data.m');
% run('./plot_data.m');



