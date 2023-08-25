#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg') # is this necessary?
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib as mpl

def discrete_sine_poly(d):
    '''
    Returns a function representing the trig-approximating polynomial of odd degree d>=3, defined for
    abs(x) <= (d+5)/2.
    
    If abs(x) > (d+5)/2, then the function representing the polynomial returns None.
    '''
    if d%2 == 0:
        raise ValueError("degree d must be odd!")
    
    def p(x):
        if abs(x) > (d+5)//2:
            return None # Outside defined range
        
        val = [0,1,1,0,-1,-1][x%6]

        if d%4 == 3: # Flip if d requires it
            val = -val
        
        if x == (d+3)/2:
            val+=1
        elif x == -(d+3)/2:
            val-=1
        elif x == (d+5)/2:
            val+=(d+2)
        elif x == -(d+5)/2:
            val-=(d+2)

        return val
    return p

def henon(p,X,x_coeff=-1):
    '''
    Given a polynomial p and X = (x,y), returns the vector (y,x_coeff*x + p(y)).

    If p(y) returns None, then this function also returns None.
    '''
    result = p(X[1])
    if result == None: # If the polynomial is not defined here, return None.
        return None
    
    return [X[1],x_coeff*X[0]+p(X[1])]

def trace_pt(p,X,escape_radius,x_coeff=-1):
    '''
    Returns the periodic orbit of a point X = (x,y) under iteration of the Henon map (x,y) -> (y,x_coeff*x + p(y)).
    Does not include the first point again.
    
    If the point is not part of a periodic cycle, returns an empty list.
    
    Parameters:
        p (function): The polynomial. Should be defined on relevant integer values.
        X (2-tuple): The starting point.
        escape_radius (int): If the point iterates outside the box of radius escape_radius, then return [].
        x_coeff (int): (Default: -1) The Henon map is (x,y) -> (y,x_coeff*x + p(y)).
    '''
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = henon(p,X,x_coeff=x_coeff)
        
        if X == None:
            return []
        if abs(X[0])>escape_radius or abs(X[1])>escape_radius:
            return []
    return orbit

def plot_orbit(orbit,colour_parameter,figure_scale,colour_style="DEFAULT"):
    '''
    Plots a given periodic cycle on the currently open figure, given some parameters.

    Parameters:
        orbit (list/tuple): The periodic cycle to plot, including each element in order, once.
                            A 2-cycle should have a list of length 2, for example.
        colour_parameter (int): A value which works together with colour_style to 
                                determine the colour to plot with. See colour_style.
        figure_scale (int): The larger this is, the smaller the circlular arrow that represents a fixed point is.
        colour_style (str): (default: "DEFAULT") Options to choose the colour of the plot.
            "DEFAULT": A collection of arbitrary colours.
            "PARAMETER": Using the cool cmap in matplotlib. 0 to 255.
    '''
    if colour_style == "DEFAULT":
        # initialise the colours
        colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
        colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
        col = colours[colour_parameter % len(colours)]         # pick the corresponding colour
    if colour_style == "PARAMETER":
        col = cm.cool(colour_parameter)

    if len(orbit) == 1:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=30/figure_scale,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit))]
        ys = [orbit[i][1] for i in range(len(orbit))]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices
        orbit.append(orbit[0])
        for k in range(len(orbit)-1):
            plt.arrow(orbit[k][0],orbit[k][1],orbit[k+1][0]-orbit[k][0],orbit[k+1][1]-orbit[k][1],width=.01,color = col,alpha =0.2)         # plots an arrow between two consecutive iterates

def find_longest_cycle_length(p,radius,x_coeff=-1):
    '''
    Returns the length of the longest periodic cycle of the Henon map (x,y) -> (y,x_coeff*x + p(y)).
    Checks only periodic cycles starting within the box of radius check_radius, and entirely contained within the box
        escape_radius.
    
    Parameters:
        p (function): The polynomial associated with the Henon map
        radius (int): The radius of the box within which we want to check for periodic cycles.
                      If the point iterates outside the box of this radius, then forget about it.
        x_coeff (int): (Default: -1) The Henon map is (x,y) -> (y,x_coeff*x + p(y)).
    '''

    found_points = []
    largest_orbit_size = 0

    for x_tocheck in range(-radius,radius):
        for y_tocheck in range(-radius,radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],radius,x_coeff=x_coeff)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
                if len(orbit) > largest_orbit_size:
                    #print(orbit)
                    largest_orbit_size = len(orbit)
    return largest_orbit_size

def create_henon_graphic(p,radius
                         ,figure_name="output"
                         ,figure_title=""
                         ,figure_size=10
                         ,reference_box_size=0
                         ,colour_style="DEFAULT"
                         ,x_coeff=-1):
    '''
    Plot all periodic cycles of the Henon map (x,y) -> (y,x_coeff*x + p(y)),
    checking the orbits of all points in a box of radius check_radius,
    and assuming all points outside of a box of radius escape_radius are not periodic (i.e. escape).
    
    Parameters:
        p (function): The polynomial associated with the Henon map
        radius (int): The radius of the box within which we want to check for periodic cycles.
        figure_name (str): (Default: "output") The name of the file to save the figure to (extension is assumed to be .png)
        figure_title (str): (Default: "") The title of the figure (shown on the plot). By default, this is empty.
        figure_size (int): (Default: 10) The size of the figure
        reference_box_size (int): (Default: 0) The size of the reference box to draw. If the value is 0, then don't draw a reference box.
        colour_style (str): (default: "DEFAULT") Options to choose the colour of the plot.
            "DEFAULT": A collection of arbitrary colours.
            "LENGTH": The colour of each cycle is determined by its length, with the longest cycle being magenta and a 1-cycle being cyan.
            "LONGEST": The colour of each cycle is determined by whether it is the longest cycle or not. 
                       Ties are handled by colouring a single cycle in pink, and the rest in purple.
        x_coeff (int): (Default: -1) The Henon map is (x,y) -> (y,x_coeff*x + p(y)).
    '''

    found_points = []        # store all the vertices whose orbits we've already plotted
    fig,ax = plt.subplots(figsize = (figure_size,figure_size))

    if reference_box_size>0: # Draw the reference box
        xs = [-reference_box_size,-reference_box_size,reference_box_size,reference_box_size,-reference_box_size]     
        ys = [-reference_box_size,reference_box_size,reference_box_size,-reference_box_size,-reference_box_size]
        plt.plot(xs,ys,"--",color = "grey")

    if colour_style == "LENGTH" or "LONGEST":
        longest_cycle = find_longest_cycle_length(p,radius,x_coeff=x_coeff)
    if colour_style == "LONGEST":
        exceptional_cycle_already_plotted = False

    count = 0
    for i in range(-radius,radius+1):    
        for j in range(-radius,radius+1):        # iterate through all the points
            if [i,j] in found_points: continue            # only iterate if we haven't plotted already, to reduce computation
            
            orbit = trace_pt(p,[i,j],radius,x_coeff=x_coeff)     # get the orbit by tracing the point
            if len(orbit) == 0: continue
            #print(orbit)
            found_points.extend(orbit)          # add each iterate to the list of plotted vertices
            if colour_style == "DEFAULT":
                plot_orbit(orbit,count,figure_scale=radius/figure_size)       # plot the orbit    
            elif colour_style == "LENGTH":
                colour_param = round((len(orbit)-1)/longest_cycle*255)
                plot_orbit(orbit,colour_param,figure_scale=radius/figure_size,colour_style="PARAMETER")
            elif colour_style == "LONGEST":
                if len(orbit) == longest_cycle:
                    if exceptional_cycle_already_plotted:
                        colour_param = 127
                    else:
                        colour_param = 255
                        exceptional_cycle_already_plotted = True
                else:
                    colour_param = 0
                plot_orbit(orbit,colour_param,figure_scale=radius/figure_size,colour_style="PARAMETER")
            count += 1                   # for count-based colour styles

    # plot formatting:
    
    if colour_style == "LENGTH":
        cmap = mpl.cm.cool
        norm = mpl.colors.Normalize(vmin=1, vmax=longest_cycle)
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax
                 ,location='bottom', label='Length of cycle', ticks=[1,longest_cycle]
                 ,fraction=0.05, pad=0.07
                 )
    elif colour_style == "LONGEST":
        cmap = mpl.cm.cool
        bounds = [1, 8, 9, 10]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
                    ,ax=ax, location='bottom', spacing='proportional'
                    ,fraction=0.05, pad=0.07)
        cbar.ax.tick_params(length=0)
        cbar.set_ticks(ticks=[1,4.5,8,9,10],labels=["","Not longest cycle","","Longest cycle(s)",""])


    plt.axis('equal')
    if figure_title != "":
        plt.title(figure_title)
    plt.savefig(figure_name)
    plt.close()

def make_your_own_function(values,index_start): 
    '''
    Define a function that takes specified values at consecutive integers starting at index_start.

    Parameters:
        values (list): The values to take at each consecutive integer starting at index_start.
        index_start (int): The first integer to define the function at.
    '''
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            #print("Outside function range!")
            return None
        return values[x-index_start]
    return p

def count_cycle_lengths(p,radius,x_coeff=-1):
    '''
    Returns a dictionary, which counts the number of periodic cycles of each cycle length, of the Henon map
    (x,y) -> (y,x_coeff*x + p(y)).

    Parameters: 
        p (function): The polynomial associated with the Henon map
        radius (int): The radius of the box within which the program will check for periodic cycles.
        x_coeff (int): (Default: -1) The Henon map is (x,y) -> (y,x_coeff*x + p(y)).
    '''
    found_points = []
    lengths = {}
    for x_tocheck in range(-radius,radius):
        for y_tocheck in range(-radius,radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],radius,x_coeff=x_coeff)
                orbit_length = len(orbit)
                if orbit_length == 0:
                    continue
                #print(orbit)
                if orbit_length not in lengths:
                    lengths[orbit_length] = 1
                else:
                    lengths[orbit_length] += 1
                
                found_points.extend(orbit)
    lengths = {a:b for a,b in sorted(lengths.items())}
    return lengths

def count_preper(p,radius,x_coeff):
    '''
    Returns the total number of integer periodic points of the Henon map (x,y) -> (y,x_coeff*x + p(y))
    within the box of the specified radius.
    
    Parameters:
        p (function): The polynomial associated with the Henon map
        radius (int): The box to check for preperiodic points.
        x_coeff (int): (Default: -1) The Henon map is (x,y) -> (y,x_coeff*x + p(y)).
    '''
    found_points = []

    for x_tocheck in range(-radius,radius):
        for y_tocheck in range(-radius,radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],radius,x_coeff=x_coeff)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
    found_points = [list(tupl) for tupl in {tuple(item) for item in found_points }]   # sanity check, ensures no duplicates     
    return len(found_points)
