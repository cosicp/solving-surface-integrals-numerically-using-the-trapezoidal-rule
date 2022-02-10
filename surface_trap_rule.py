# Python 3 code for solving surface integrals numerically
# author: Petar Cosic
# email: cpetar112@gmail.com
# github: https://github.com/cosicp
# date: 10.2.2022.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def trap2D(x, yg, yd, F):
    """function that solves and plots the solution of a surface integral for a given domain"""
    nx = int(len(x)-1)
    ny = int((len(x)-1)/2)
    xtick = nx+1
    ytick = ny+1
    a = x[0]
    b = x[xtick-1]
    hy = (yg - yd)/ ny
    hx = (b - a) / nx  
    
    # sovling the edges of the surface
    Part1 = 1/4*hx*hy[0]*F(a,yd[0])
    Part2 = 1/4*hx*hy[0]*F(a,yg[0])
    Part4 = 1/4*hx*hy[xtick-1]*F(b,yd[xtick-1])
    Part5 = 1/4*hx*hy[xtick-1]*F(b,yg[xtick-1])

    # solving line 1/4
    sum_3 = 0
    yi_1=yd[0]
    yi_end=yd[xtick-1]
    for i in range(1,ny):
        sum_3=sum_3+F(a,yi_1+(i)*hy[0])
    Part3 = 1/4*hx*hy[0]*2*sum_3

    # solving line 2/4
    sum_6 = 0
    yi_1=yd[0]
    yi_end=yd[xtick-1]
    for i in range(2,ny):
        sum_6=sum_6+F(b,yi_end+(i)*hy[xtick-1])
    Part6 = 1/4*hx*hy[xtick-1]*2*sum_6
    
    # solving line 3/4
    sum_7 = 0
    for j in range(1,nx):
        sum_7 = sum_7 + 1/2*hy[j]*F(x[j],yd[j])
    Part7 = hx*sum_7
    
    # solving line 4/4
    sum_8 = 0
    for j in range(1,nx):
        sum_8=sum_8 + (1/2)*hy[j]*F(x[j],yg[j])
    Part8 = hx*sum_8    

    # solving the inner part 
    sum_9_1 = 0
    for j in range(1,nx):
        sum_9_2 = 0
        for i in range(1, ny):
            sum_9_2 = sum_9_2 + F(x[j], (yd[j] + (i)*hy[j])) 
        sum_9_1 = sum_9_1 + sum_9_2*hy[j]
    Part9 = sum_9_1 * hx

    # plotting the solution
    Y = 0*np.ones([xtick,ytick])
    for j in range(0,xtick):
        for i in range(0,ytick):
            Y[j,i]=yd[j]+i*hy[j]

    X = np.ones([xtick, ytick])
    for j in range(0,xtick):
        X[j,:]=X[j,:]*x[j]
        
    X = np.transpose(X)
    Y = np.ones([ytick,xtick])
    for i in range(0,ytick):
        for j in range(0,xtick):
            Y[i,j] = Y[i,j]*(yd[j]+(i)*hy[j])
   
    Z1= np.sqrt(100-X**2-Y**2+0j)
    Z = np.real(Z1)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(X, Y, Z ,cmap=plt.cm.cividis,linewidth=0, antialiased=False)
    ax.grid(False)
    ax.set_xlabel('$x$', fontsize=12)
    ax.set_ylabel('y',fontsize=12)
    ax.set_zlabel('$z = \sqrt{100-x^2-y^2}$',fontsize=10)
    ax.set_title(f'num of nodes ={len(x)}',fontsize=10)
    for ii in range(0,360,30):
       ax.view_init(elev=10., azim=ii)
       plt.savefig("images/img%d.svg" % ii)
    plt.show()
    
    return Part1+Part2+Part3+Part4+Part5+Part6+Part7+Part8+Part9

x = np.linspace(0,10,101)
yd = 0*x
yg = np.sqrt(100-x**2)
F = lambda x,y: np.sqrt(100-x**2-y**2+0j)
sol  = trap2D(x, yg, yd, F)
print(f'The numerical solution for the integral of the function sqrt(100-x**2-y**2)\
    over the domain x in (0,10) , y in (0,sqrt(100-x^2)) is {np.real(sol):.3f}')

