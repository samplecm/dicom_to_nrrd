import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math
from Contours import Contours
from numpy.lib.type_check import imag
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider

from PIL import Image
import copy

import plotly.graph_objects as go

def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.
    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

def plotSubsegments(structure):
    print("In PlotStructure")
    fig = plt.figure()
    ax : plt.Axes = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)


    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    minX = 1000
    minY = 1000
    minZ = 1000
    maxX = -1000
    maxY = -1000
    maxZ = -1000


    colours = ['r', 'b', 'g', 'y', 'm', 'c', 'k']


    for i in range(len(structure)):
        substructure = structure[i]
        # colour = colours[colour_idx]
        # colour_idx = (colour_idx + 1) % 7
        colour = colours[i % 7]
        for c, contour in enumerate(substructure):   
            if len(contour) == 0: 
                continue  
            if len(contour[0]) == 0:
                continue    
            # if c % 2 != 0:
            #     continue
            x = []
            y = []
            z = []
            for point in contour[0]:
                if point[0] > maxX:
                    maxX = point[0]
                elif point[0] < minX:
                    minX = point[0]
                if point[1] > maxY:
                    maxY = point[1]
                if point[1] < minY:
                    minY = point[1]
                if point[2] > maxZ:
                    maxZ = point[2]
                if point[2] < minZ:
                    minZ = point[2]                     
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])

            #ax.plot(x,y,z, colour)  
            points = [list(zip(x,y,z))]
            poly = a3.art3d.Poly3DCollection(points)  
            poly.set_color(colour)
            poly.set_edgecolor('k')
            ax.add_collection3d(poly)

    

    ax.set_xlim((minX-5, maxX+5))    
    ax.set_ylim((minY-5, maxY+5))    
    ax.set_zlim((minZ-5, maxZ+5))    

    set_axes_equal(ax)
    ax.grid(False) 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([]) 
    plt.axis('off')
    plt.show()
    print("")

def Get_Colormap_RGB(val):
    
    # val_blue = min(0.5,val)
    # blue = (1-2*val_blue)

    # val_green = abs(val-0.5)
    # green = (1-2*val_green)

    # val_red = 1 - max(0.5,val)
    # red = (1-2*val_red)
    # i reveresed the order by accident... so reverse val
    val = 1 - val

    if val < 0.25:
        red = 1
        blue = 0.1
        green = 0.9*(val/0.25) + 0.1
    elif val < 0.5:
        red = -0.9*((val-0.25)/0.25) + 1
        green = 1
        blue = 0.1
    elif val < 0.75:
        red = 0.1
        green = 1    
        blue = 0.9*((val-0.5)/0.25) + 0.1
    else:
        red = 0.1
        blue = 1
        green = -0.9*((val-0.75)/0.25) + 1


    return (red, green, blue, 1)
    
def plotStructure_unfilled(structure):
    print("In PlotStructure")
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    minX = 1000
    minY = 1000
    minZ = 1000
    maxX = -1000
    maxY = -1000
    maxZ = -1000

    colour_idx = 0
    colours = ['r', 'b', 'g', 'y', 'm', 'c', 'k']
    for i in range(len(structure)):
        substructure = structure[i]
        colour = colours[colour_idx]
        colour_idx = (colour_idx + 1) % 7
        for contour in substructure:          
            x = []
            y = []
            z = []
            for point in contour:
                if point[0] > maxX:
                    maxX = point[0]
                elif point[0] < minX:
                    minX = point[0]
                if point[1] > maxY:
                    maxY = point[1]
                if point[1] < minY:
                    minY = point[1]
                if point[2] > maxZ:
                    maxZ = point[2]
                if point[2] < minZ:
                    minZ = point[2]                     
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])
            ax.plot(x,y,z, colour)    
    ax.set_xlim((minX-5, maxX+5))    
    ax.set_ylim((minY-5, maxY+5))    
    ax.set_zlim((minZ-5, maxZ+5))    

        
    plt.show()
    print("")