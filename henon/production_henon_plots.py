#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib as mpl

def discrete_sine_poly(d,negative = False):      # returns the discrete sine polynomial p.
                                # defined in a box, outside of which the Henon maps diverge to infty.
    if d%4==1: # given by 0 1 1 0 -1 -1 0 ... (starting at 0)
        def p(x):
            # define it periodically 
            if x%3 == 0:
                val = 0
            elif x%6 == 1 or x%6 == 2:
                val = 1
            else:
                val = -1

            # adjust the endpoints and return the appropriate value of p(x)
            if abs(x)<=(d+1)/2:
                pass
            elif x == (d+3)/2:
                val+=1
            elif x == -(d+3)/2:
                val-=1
            elif x == (d+5)/2:
                val+=(d+2)
            elif x == -(d+5)/2:
                val-=(d+2)
            else:  # can only calculate for x in the "nice" range
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
            
            if negative:
                val = -val
            return val
        return p
    elif d%4==3: # given by 0 -1 -1 0 1 1 0 ... (starting at 0)
        def p(x):
            # define it periodically 
            if x%3 == 0:
                val = 0
            elif x%6 == 1 or x%6 == 2:
                val = -1
            else:
                val = 1
                
            # adjust the endpoints and return the appropriate value of p(x)
            if abs(x)<=(d+1)/2:
                pass
            elif x == (d+3)/2:
                val+=1
            elif x == -(d+3)/2:
                val-=1
            elif x == (d+5)/2:
                val+=(d+2)
            elif x == -(d+5)/2:
                val-=(d+2)
            else: # can only calculate for x in the "nice" range
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
            
            if negative:
                val = -val
            return val
        return p
    
    print("Warning, this value of d has not yet been implemented!")
    return None

def henon(p,X,x_coefficient=-1):     # this is the Henon map of polynomial p: (x,y) -> (y, x_coefficient*x + p(y))
    result = p(X[1])
    if result == None:
        # print("Henon map failed at X=",X)
        return None # This is a (carried over) flag telling us that we went outside of the nice box, and hence we should stop iterating
    else:
        return [X[1],x_coefficient*X[0]+p(X[1])]

def trace_pt(p,X,box_range,x_coefficient=-1):        # follows the orbit of a point X under Henon map with polynomial p, 
                                    # stopping either upon repeating a vertex or when we go outside the nice range
                                    # Returns the path that the point took, only if it's a loop.
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = henon(p,X,x_coefficient=x_coefficient)
        if X == None:
            return []
        if abs(X[0])>box_range or abs(X[1])>box_range:
            return []               # if we iterate outside the d+2 box, we don't want to plot the orbit at all
    return orbit

def plot_orbit (orbit,colour_parameter,box_range,colour_style="DEFAULT",figure_size=10):         # plot the orbit
    
    if colour_style == "DEFAULT":
        # initialise the colours
        colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
        colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
        col = colours[colour_parameter % len(colours)]         # pick the corresponding colour
    if colour_style == "PARAMETER":
        col = cm.cool(colour_parameter)

    if len(orbit) == 1:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=30*figure_size/box_range,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit))]
        ys = [orbit[i][1] for i in range(len(orbit))]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices
        orbit.append(orbit[0])
        for k in range(len(orbit)-1):
            plt.arrow(orbit[k][0],orbit[k][1],orbit[k+1][0]-orbit[k][0],orbit[k+1][1]-orbit[k][1],width=.01,color = col,alpha =0.2)         # plots an arrow between two consecutive iterates

def find_longest_cycle_length(p,escape_radius,check_radius,x_coefficient=-1): # find the length of the longest cycle in an orbit p.
                                                             # can easily be modified to give you the actual cycle too.
    found_points = []
    largest_orbit_size = 0

    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius,x_coefficient=x_coefficient)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
                if len(orbit) > largest_orbit_size:
                    #print(orbit)
                    largest_orbit_size = len(orbit)
    return largest_orbit_size

def create_henon_graphic(p,escape_radius,check_radius,
                         figure_name="output",figure_title="",
                         figure_size=10,
                         reference_box_size=0,
                         colour_style="DEFAULT",
                         x_coefficient=-1):

    
    found_points = []        # store all the vertices whose orbits we've already plotted
    fig,ax = plt.subplots(figsize = (figure_size,figure_size))

    if reference_box_size>0: # Draw the reference box
        xs = [-reference_box_size,-reference_box_size,reference_box_size,reference_box_size,-reference_box_size]     
        ys = [-reference_box_size,reference_box_size,reference_box_size,-reference_box_size,-reference_box_size]
        plt.plot(xs,ys,"--",color = "grey")
    # plt.grid(which="both")
    # plt.xticks([i for i in range(-box_range-1,box_range+1)])
    # plt.yticks([i for i in range(-box_range-1,box_range+1)])

    if colour_style == "LENGTH" or "LONGEST":
        longest_cycle = find_longest_cycle_length(p,escape_radius,check_radius,x_coefficient=x_coefficient)
    if colour_style == "LONGEST":
        exceptional_cycle_already_plotted = False

    count = 0
    for i in range(-check_radius,check_radius+1):    
        for j in range(-check_radius,check_radius+1):        # iterate through all the points
            if [i,j] in found_points: continue            # only iterate if we haven't plotted already, to reduce computation
            
            orbit = trace_pt(p,[i,j],escape_radius,x_coefficient=x_coefficient)     # get the orbit by tracing the point
            if len(orbit) == 0: continue
            #print(orbit)
            found_points.extend(orbit)          # add each iterate to the list of plotted vertices
            if colour_style == "DEFAULT":
                plot_orbit(orbit,count,escape_radius,figure_size=figure_size)       # plot the orbit    
            elif colour_style == "LENGTH":
                colour_param = round((len(orbit)-1)/longest_cycle*255)
                plot_orbit(orbit,colour_param,escape_radius,colour_style="PARAMETER",figure_size=figure_size)
            elif colour_style == "LONGEST":
                if len(orbit)-1 == longest_cycle:
                    if exceptional_cycle_already_plotted:
                        colour_param = 127
                    else:
                        colour_param = 255
                        exceptional_cycle_already_plotted = True
                else:
                    colour_param = 0
                plot_orbit(orbit,colour_param,escape_radius,colour_style="PARAMETER",figure_size=figure_size)
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
    ax.set_xticks(ticks=list(range(-100,101)))
    ax.set_axisbelow(True)
    plt.grid(which="both",axis="both",color='#CCCCCC')
    plt.savefig(figure_name)
    plt.close()

def make_your_own_function(values,index_start): # define a function that takes specified values at integer intervals starting at index_start.
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            #print("Outside function range!")
            return None
        return values[x-index_start]
    return p

poly = discrete_sine_poly(3)
create_henon_graphic(poly,6,6,colour_style="DEFAULT",figure_size=5)