import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

def plot_1():
    def f1(x):
        return 1/10*x**3 + 1/20*x**2 - 3/10*x
    def f2(x):
        # return ((x+4)**2)*(x<0) + ((x-4)**2)*(x>=0)
        return f1(x+4)*(x<0) + f1(-x+4)*(x>=0)
    x1 = np.linspace(-2, 2, 5)
    x2 = np.linspace(-2, 2, 5)
    x3 = np.append(np.linspace(-6, -2, 5), np.linspace(2, 6, 5))
    xtilde1 = 1.25
    xtilde2 = 3.5
    xtilde3 = 0.0
    
    y1 = f1(x1)
    y2 = f1(x2)
    y3 = f2(x3)
        
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 3.5), sharex=True)
    ax1.plot(x1, y1, marker='o', linewidth=0)
    ax1.plot([xtilde1, xtilde1], [-10, 10], c='r')
    ax1.set_yticks([])
    ax1.set_ylim([-0.7, 1])
    ax2.plot(x2, y2, marker='o', linewidth=0)
    ax2.plot([xtilde2, xtilde2], [-10, 10], c='r')
    ax2.set_yticks([])
    ax2.set_ylim([-0.7, 1])
    ax3.plot(x3, y3, marker='o', linewidth=0)
    ax3.plot([xtilde3, xtilde3], [-10, 10], c='r')
    ax3.set_yticks([])
    ax3.set_ylim([-0.7, 1.0])
    fig.tight_layout()


def plot_2():
    np.random.seed(13)
    points = np.random.normal(size=(30,2))
    hull = ConvexHull(points)
    
    hull_vertices = hull.vertices
    hull_points = points[hull_vertices]
    
    plt.figure(figsize=(4, 4))
    plt.fill(hull_points[:, 0], hull_points[:, 1], 'lightblue', alpha=0.5)
    plt.scatter(points[:, 0], points[:, 1], label='points')
    # plt.scatter(hull_points[:, 0], hull_points[:, 1], color='red', label='hull vertices')
    # plt.legend()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    plt.tight_layout()

if __name__ == '__main__':
    plot_1()
    plot_2()