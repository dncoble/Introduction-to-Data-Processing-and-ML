import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
"""
FEA representation of transverse motion of a fixed-free Euler-Bernoulli beam.
I use the naming conventions of An Introduction to the Finite Element
Methodology by JN Reddy. 
"""
def fea_beam(f, t=200):
    #%% description of the beam
    n_nodes = 40
    EI = 1
    rhoA = 2
    L = 1
    node_points = [L*i/(n_nodes-1) for i in range(n_nodes)]

    deltat = .01
    t = 200 # total time
    #%% related values
    n_elements = n_nodes - 1
    t_timesteps = int(t/deltat)
    #%% create finite element matrices K, M
    # create global equations
    K = np.zeros((2*n_nodes,2*n_nodes))
    M = np.zeros((2*n_nodes,2*n_nodes))
    for i in range(1, n_elements+1):
        h = node_points[i] - node_points[i-1]
        K_el = 2*EI/h**3*\
            np.matrix([[6, -3*h, -6, -3*h],
                      [-3*h, 2*h**2, 3*h, h**2],
                      [-6, 3*h, 6, 3*h],
                      [-3*h, h**2, 3*h, 2*h**2]])
        M_el = rhoA*h/420*\
            np.matrix([[156, -22*h, 54, 13*h],
                      [-22*h, 4*h**2, -13*h, -3*h**2],
                      [54, -13*h, 156, 22*h],
                      [13*h, -3*h**2, 22*h, 4*h**2]])
        
        K[2*i-2:2*i+2,2*i-2:2*i+2] += K_el
        M[2*i-2:2*i+2,2*i-2:2*i+2] += M_el
    Minv = np.linalg.inv(M)
    #%% do finite element, discrete time
    """
    Newmark's Scheme
    alpha = 3/2, gamma = 8/5 is the Galerkin method
    """
    alpha = 3/2
    gamma = 8/5
    # related values
    beta = gamma/2
    a1 = alpha*deltat
    a2 = (1 - alpha)*deltat
    a3 = 1/(beta*deltat**2)
    a4 = a3*deltat
    a5 = 1/gamma - 1
    # a6 = alpha/(beta*deltat)
    # a7 = alpha/beta - 1
    # a8 = (alpha/gamma - 1)*deltat
    Khat = K + a3*M
    Khatsub = Khat[2:,2:]
    Khatsubinv = np.linalg.inv(Khatsub)
    # displacement over time [timesteps, node]
    w_t = np.zeros([t_timesteps, 2*n_nodes])
    # initial values
    F_s = np.zeros(2*n_nodes); F_s[2*n_nodes-2] = f(0)
    w_s = np.zeros(2*n_nodes); w_t[0] = w_s
    dotw_s = np.zeros(2*n_nodes)
    ddotw_s = Minv@(F_s - K @ w_s)
    for s in range(1, t_timesteps):
        t = s*deltat
        F_s = np.zeros(2*n_nodes); F_s[2*n_nodes-2] = f(t)
        # computation order for Newmark's Scheme
        b_s = w_s + deltat*dotw_s + .5*(1-gamma)*(deltat**2)*ddotw_s
        # B_s = a6*w_s + a7*dotw_s + a8*ddotw_s
        A_s = a3*b_s
        Fhat = F_s + M @ A_s
        # new timesteps
        w_s1 = np.append(np.zeros((2,)), Khatsubinv@Fhat[2:])
        ddotw_s1 = a3*(w_s1 - w_s) - a4*dotw_s - a5*ddotw_s
        dotw_s1 = dotw_s + a2*ddotw_s + a1*ddotw_s1
        
        w_s = w_s1
        dotw_s = dotw_s1
        ddotw_s = ddotw_s1
        
        w_t[s] = w_s
    
    t = np.arange(0, t_timesteps)*deltat
    x = np.array(node_points)
    W = w_t[:,::2]
    return (W, t, x)