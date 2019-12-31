function action = csecond_order_fun(t, x, oracle)
    persistent s_n s_s t_s;
    
    if(isempty(s_n))
        s_n = oracle(x);
    end
    
    if(isempty(s_s))
        s_s = oracle(x);
    end
    
    if(isempty(t_s))
        t_s = t;
    end
        
    s_n = s_n + oracle(x);
    s_s = [s_s; s_n];
    t_s = [t_s; t];
    
    assigin('base', 'ts', t_s);
    assigin('base', 'actions', s_s);

    action = s_n;
end