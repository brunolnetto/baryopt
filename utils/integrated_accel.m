function z_bar = integrated_accel(m_1, xhat_1, delta_xhat_1, ...
                                  lambda_z, nu, oracle)
    persistent m_2 z_bar_1;
    
    if(isempty(m_2))
        m_2 = 0;
    end
    
    if(isempty(z_bar_1))
        z_bar_1 = zeros(size(delta_xhat_1));
    end
    
    % Calculation of mean for stochastic function
    e_z = exp(-nu*oracle(xhat_1));
    z_bar = -(1/m_1)*(lambda_z*e_z*delta_xhat_1 + z_bar_1*m_2);
    
    m_2 = m_1;
    z_bar_1 = z_bar;
end