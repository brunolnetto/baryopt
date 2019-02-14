oracle = @(x) sum(x.^2);

a = -10;
b = 10;
n = 100;

x0 = [0, 0];
norms = [];
nus = 0.01:0.001:5;

for nu = nus 
    % Batch version
    x_ = uniform(a, b, [n, 2]);
    x = expbary(oracle, x_, nu);
    norms = [norms; norm(x - x0)];
end

figure
plot(nus, norms)