function action = dsecond_order_fun(x, oracle)
    persistent s_n s_s;
    
    if(isempty(s_n))
        s_n = oracle(x);
    end

    if(isempty(s_s))
        s_s = oracle(x);
    end
        
    s_n = s_n + oracle(x);
    s_s = [s_s; s_n];
    
    assignin('base', 'actions', s_s);
    action = s_n;
end