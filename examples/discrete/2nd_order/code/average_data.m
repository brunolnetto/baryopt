x_star = zeros(1, 2);
for i = 1:length(xhat_tests)
    x_aux = xstar_tests{i};
    x_star = x_star + x_aux;
end

x_star = x_star/n_iterations;

xhat_mean = zeros(size(xhat_tests{1}));
x_mean = zeros(size(xs_tests{1}));
y_mean = zeros(size(xs_tests{1}));
for i = 1:length(xhat_tests)
    x_mean = x_mean + xs_tests{i};
    xhat_mean = xhat_mean + xhat_tests{i};
    y_mean = y_mean + ys_tests{i};
end

xhat_mean = xhat_mean/n_iterations;
x_mean = x_mean/n_iterations;
y_mean = y_mean/n_iterations;

zbar_mean = zeros(size(zbars_tests{1}));
for i = 1:length(zbars_tests)
    zbar_mean = zbar_mean + zbars_tests{i};
end

zbar_mean = zbar_mean/n_iterations;

fs_mean = [];
for i = 1:iterations
    f_i = oracle(xhat_mean(i, :));
    fs_mean = [fs_mean; f_i];
end

delta_mean = zeros(size(delta_tests{1}));
for i = 1:length(delta_tests)
    delta_mean = delta_mean + delta_tests{i};
end

delta_mean = delta_mean/n_iterations;

