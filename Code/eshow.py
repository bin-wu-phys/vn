#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:23:58 2019

@author: biwu
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def XYZ(rLIST, thetaLIST, eLIST):
    thLIST = np.concatenate((thetaLIST, np.array([2.0*np.pi])))

    # Create the mesh in polar coordinates and compute corresponding Z.
    R, Th = np.meshgrid(rLIST, thLIST)
    f = np.transpose(eLIST)
    Z = np.concatenate((f, np.reshape(f[0, :], (1, len(f[0,:])))))

    # Express the mesh in the cartesian system.
    X, Y = R*np.cos(Th), R*np.sin(Th)
    return X, Y, Z

def energyContour(rLIST, thetaLIST, eLIST, t):
    X, Y, Z = XYZ(rLIST, thetaLIST, eLIST)
    
    # Canvas and plot
    fig = plt.figure(figsize=(5.0,5.0))
    ax = fig.add_subplot(111#, projection='3d'
                        )

    # Plot the surface.
    ax.contour(X, Y, Z)


    # Tweak the limits and add latex math labels.
    ax.set_title('t = {:.2f}'.format(t))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # Set x and y ranges
    upper = rLIST[-1]
    ax.set_xlim(-upper, upper)
    ax.set_ylim(-upper, upper)

    plt.grid()
    plt.show()
    return fig

def energy3D(rLIST, thetaLIST, eLIST, t):
    X, Y, Z = XYZ(rLIST, thetaLIST, eLIST)
    
    # Canvas and plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d'
                        )
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm)
    # ax.view_init(elev=70, azim=45)

    # Tweak the limits and add latex math labels.
    ax.set_title('t = {:.2f}'.format(t))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('$\epsilon$')
    #ax.set_zticks([0, 0.007, 0.014])

# Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.6, aspect=10)

    plt.show()
    return fig