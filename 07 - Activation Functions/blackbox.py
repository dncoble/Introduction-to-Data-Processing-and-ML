import numpy as np
"""
No peaking!
"""
def f(x):
    return np.sin(.5*x)+np.sin(x)+np.cos(x**2/8)+x/10

# import matplotlib.pyplot as plt
# x = np.linspace(-10, 10, num=200)
# plt.figure()
# plt.plot(x, f(x))