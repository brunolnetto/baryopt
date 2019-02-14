function xhat = expbary(oracle, x, nu)
% Batch barycenter algorithm for direct optimization
% https://arxiv.org/abs/1801.10533
% In:
%   - oracle [function]: Oracle function
%   - x []: Query values
%   - nu []: positive value (Caution on its value due overflow)
% Out:
%   - xhat []: Optimum position

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