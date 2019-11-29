function z = curiosity_fun(t, xhat, zeta)
    xhat_perp = null(xhat');
    z = zeta*xhat_perp(:, 1);
end