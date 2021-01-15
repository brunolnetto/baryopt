% Take the mean value
    
% X star
xstar = zeros(1, 2);
for i = 1:length(xhat_tests)
    x_aux = xstar_tests{i};
    xstar = x_star + x_aux;
end

xstar_mean = xstar/n_iterations;

% X hat and x
xhat_mean = zeros(size(xhat_tests{1}));
x_mean = zeros(size(xs_tests{1}));
y_mean = zeros(size(ys_tests{1}));
for i = 1:length(xhat_tests)
    x_mean = x_mean + xs_tests{i};
    xhat_mean = xhat_mean + xhat_tests{i};
    y_mean = y_mean + ys_tests{i};
end

xhat_mean = xhat_mean/n_iterations;
x_mean = x_mean/n_iterations;
y_mean = y_mean/n_iterations;

% m
m_mean = zeros(size(ms_tests{1}));
for i = 1:length(xhat_tests)
    m_mean = m_mean + ms_tests{i};
    m_mean = m_mean + ms_tests{i};
end

m_mean = m_mean/n_iterations;

% Deltas
delta_mean = zeros(size(delta_tests{1}));
for i = 1:length(zbars_tests)
    delta_mean = delta_mean + delta_tests{i};
end

delta_mean = delta_mean/n_iterations;

% zbars
zbar_mean = zeros(size(zbars_tests{1}));
for i = 1:length(zbars_tests)
    zbar_mean = zbar_mean + zbars_tests{i};
end

zbar_mean = zbar_mean/length(zbars_tests);

% function value
fs_mean = [];
for i = 1:iterations
    f_i = oracle(xhat_mean(i, :));
    fs_mean = [fs_mean; f_i];
end

xstar_means{end+1} = xstar_mean;
x_means{end+1} = x_mean;
y_means{end+1} = y_mean;
xhat_means{end+1} = xhat_mean;
zbars_means{end+1} = zbar_mean;
ms_means{end+1} = m_mean;
fs_means{end+1} = fs_mean;
deltas_means{end+1} = delta_mean;