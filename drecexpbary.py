def drecexpbary(oracle, x0,
                nu = 1, sigma = 0.1, 
                zeta = 0, lambda_ = 1, 
                iterations = 100):
    '''
     Recursive barycenter algorithm for direct optimization
     In:
       - oracle [function]: Oracle function e.g. lambda(x) numpy.power(x, 2)
       - x0 [np.array]: Initial query values
       - nu [double]: positive value (Caution on its value due overflow)
       - sigma [double]: Std deviation of normal distribution

        - zeta [double]: Proportional value for mean of normal distribution
        - lambda [double]: Forgetting factor between 0 and 1
        - iterations [int]: Maximum number of iterations
     Out:
        - x [np.array]: Optimum position
    '''
    import numpy as np
    import scipy as scp

    # Initialization
    xhat_1 = x0
    m_1 = 0

    len_x0 = x0.__len__()
    deltax_1 = np.zeros((len_x0, 1))
    solution_found = False

    normrnd = scp.random.normal

    # Optimization loop
    i = 1
    while(not solution_found):
        i = i + 1

        z = normrnd(zeta*deltax_1, sigma).T
        
        x = xhat_1 + z
        fi = oracle(x)

        e_i = np.exp(-nu*fi)
        m = lambda_*m_1 + e_i
        xhat = (1/m)*(lambda_*m_1*xhat_1 + x*e_i)
        
        solution_found = i >= iterations
        
    return xhat
