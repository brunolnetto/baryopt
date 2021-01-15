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
x_0 = 0;
y_0 = 0;
oracle = @(x) (x(1) - x_0)^2 + (x(2) + y_0)^2;

n_points = 10;
n_batches = 10;

batches = {};
xhats = {};
Ss = {};

for i = 1:n_batches
    size_batch = [n_points, 2];
    batch_i = rand_ab(-1, 1, [n_points, 2]);
    
    batches{end+1} = batch_i;
    
    [xhat, Si] = expbary(oracle, batch_i, nu);
    
    xhats{end+1} = xhat;
    Ss{end+1} = Si;
end

acc = zeros(1, 2);
S_final = 0;
for i = 1:n_batches
    acc = acc + Ss{i}*xhats{i};
    S_final = S_final + Ss{i};
end

xhat_final = acc/S_final;

FontSize = 20;

hfig = my_figure();

for i = 1:n_batches
    plot(batches{i}(:, 1), ...
    batches{i}(:, 2), '*');
    hold on;
    xhat_i = xhats{i};
    
    plot(xhat_i(1),xhat_i(2) , '.', 'MarkerSize', 30);
end

hold on

plot(xhat_final(1), xhat_final(2) , ...
     'kD', 'MarkerSize', 12, ...
     'MarkerFaceColor','green');

hold on

plot(x_0, y_0 , ...
     '-p', 'MarkerSize', 12, ...
     'MarkerFaceColor','red'); 
 
axis square

title('Evaluation of the barycenter method for multiple batches', ...
      'interpreter', 'latex', 'FontSize', FontSize);

set(gca,'FontSize', FontSize);

% Save folder
path = [pwd '/../imgs/'];
fname = 'multiple_points';
saveas(hfig, [path, fname], 'epsc')

