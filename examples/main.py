import numpy as np

if __name__ == "__main__":
	# Two dimensional y = x1^2 + x2^2
    oracle = lambda(x): np.power(x, 2)
    
    # Initial conditions
    x0 = np.array([0, 0])
    
    # Function call
    x = drecexpbary(oracle, x0, 10, 0.1, 0, 1, 100)
    print(x)
