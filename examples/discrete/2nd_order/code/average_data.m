x_star = zeros(1, 2);
for i = 1:length(xhat_tests)
    x_aux = xstar_tests{i};
    x_star = x_star + x_aux;
end

x_star = x_star/n_iterations;

xhat_mean = zeros(size(xhat_tests{1}));
vhat_mean = zeros(size(vhat_tests{1}));
x_mean = zeros(size(xs_tests{1}));

for i = 1:length(xhat_tests)
    x_mean = x_mean + xs_tests{i};
    xhat_mean = xhat_mean + xhat_tests{i};
    vhat_mean = vhat_mean + vhat_tests{i};
end

vhat_mean = vhat_mean/n_iterations;
xhat_mean = xhat_mean/n_iterations;
x_mean = x_mean/n_iterations;

zbar_mean = zeros(size(zbars_tests{1}));
for i = 1:length(zbars_tests)
    zbar_mean = zbar_mean + zbars_tests{i};
end

zbar_mean = zbar_mean/n_iterations;

fs_mean = [];
for i = 1:iterations
    vhat = vhat_mean(i, :);
    f_i = oracle(xhat_mean(i, :)) + vhat*vhat';
    fs_mean = [fs_mean; f_i];
end

