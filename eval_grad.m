syms x y
jac_func = jacobian(func, [x; y]);

[m, ~] = size(xs);

wb = my_waitbar('Loading gradients...');
grad_vals = [];
for i = 1:m
    x_i = xs(i, :);
    grad_val = double(subs(jac_func, [x, y], x_i));
    grad_vals = [grad_vals; grad_val/norm(grad_val)];

    wb.update_waitbar(i, m);
end

wb = my_waitbar('Loading Deltas...');
delta_xhat_ns = [];
for i = 2:m
    delta_x_n = xs(i, :) - xs(i - 1, :);
    delta_xhat_ns = [delta_xhat_ns; delta_x_n/norm(delta_x_n)];

    wb.update_waitbar(i, m);
end

wb = my_waitbar('Loading projections...');
dots_delta_grad = [];
for i = 1:m-1
    delta_xhat_n = delta_xhat_ns(i, :);
    grad_val = grad_vals(i, :);

    dot_dg = dot(delta_xhat_n, grad_val);
    dots_delta_grad = [dots_delta_grad; dot_dg];

    wb.update_waitbar(i, m-1);
end

histogram(dots_delta_grad, 100);
