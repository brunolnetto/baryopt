function xhat = expbary(oracle, x, nu)
    [m, n] = size(x);
    
    Si = 0;
    for i = 1:m
        e_ = exp(-nu*oracle(x(i, :)));
        Si = Si + e_;
    end
    
    Sxi = zeros(1, n);
    for i = 1:m
        Sxi = Sxi + x(i, :)*exp(-nu*oracle(x(i, :)));
    end
    
    xhat = Sxi/Si;
end