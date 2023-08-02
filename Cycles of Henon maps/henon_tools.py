#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

def discrete_sine_poly(d):      # returns the discrete sine polynomial p.
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
                return val
            elif x == (d+3)/2:
                return val+1
            elif x == -(d+3)/2:
                return val-1
            elif x == (d+5)/2:
                return val + (d+2)
            elif x == -(d+5)/2:
                return val - (d+2)
            else:  # can only calculate for x in the "nice" range
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
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
                return val
            elif x == (d+3)/2:
                return val+1
            elif x == -(d+3)/2:
                return val-1
            elif x == (d+5)/2:
                return val + (d+2)
            elif x == -(d+5)/2:
                return val - (d+2)
            else: # can only calculate for x in the "nice" range
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
    else:
        print("Warning, this value of d has not yet been implemented!")
        return None
    return p
            
def henon(p,X):     # this is the Henon map of polynomial p: (x,y) -> (y, -x + p(y))
    result = p(X[1])
    if result == None:
        # print("Henon map failed at X=",X)
        return None # This is a (carried over) flag telling us that we went outside of the nice box, and hence we should stop iterating
    else:
        return [X[1],-X[0]+p(X[1])]

def trace_pt(p,X,box_range):        # follows the orbit of a point X under Henon map with polynomial p, 
                                    # stopping either upon repeating a vertex or when we go outside the nice range
                                    # Returns the path that the point took, only if it's a loop.
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = henon(p,X)
        if abs(X[0])>box_range or abs(X[1])>box_range:
            return []               # if we iterate outside the d+2 box, we don't want to plot the orbit at all
    orbit.append(X)                 # add the last vertex to complete the orbit
    return orbit

def plot_orbit (orbit,colour_index,box_range):         # plot the orbit
    # initialise the colours
    colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
    colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
    col = colours[colour_index % len(colours)]         # pick the corresponding colour


    if len(orbit) == 2:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=300/box_range,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit)-1)]
        ys = [orbit[i][1] for i in range(len(orbit)-1)]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices

        for k in range(len(orbit)-1):
            plt.arrow(orbit[k][0],orbit[k][1],orbit[k+1][0]-orbit[k][0],orbit[k+1][1]-orbit[k][1],width=.01,color = col,alpha =0.2)         # plots an arrow between two consecutive iterates

def find_largest_orbit(p,escape_radius,check_radius):
    found_points = []
    largest_orbit_size = 0

    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
                if len(orbit)-1 > largest_orbit_size:
                    #print(orbit)
                    largest_orbit_size = len(orbit)-1
    return largest_orbit_size

def create_henon_graphic(p,escape_radius,check_radius,figure_name="output",reference_box_size=0):
    plotted = []        # store all the vertices whose orbits we've already plotted

    colour_index = 0
    plt.figure(figsize = (15,15))

    if reference_box_size>0: # Draw the reference box
        xs = [-reference_box_size,-reference_box_size,reference_box_size,reference_box_size,-reference_box_size]     
        ys = [-reference_box_size,reference_box_size,reference_box_size,-reference_box_size,-reference_box_size]
        plt.plot(xs,ys,"--",color = "grey")
    # plt.grid(which="both")
    # plt.xticks([i for i in range(-box_range-1,box_range+1)])
    # plt.yticks([i for i in range(-box_range-1,box_range+1)])

    for i in range(-check_radius,check_radius+1):    
        for j in range(-check_radius,check_radius+1):        # iterate through all the points
            if [i,j] not in plotted:            # only iterate if we haven't plotted already, to reduce computation
                orbit = trace_pt(p,[i,j],escape_radius)     # get the orbit by tracing the point
                for pt in orbit:
                    plotted.append(pt)          # add each iterate to the list of plotted vertices
                plot_orbit(orbit,colour_index,escape_radius)       # plot the orbit    
                colour_index +=1                   # to vary our colours in a nice way

    # plot formatting:

    plt.axis('equal')
    plt.savefig(figure_name)
    plt.close()

def print_largest_orbits_discrete_sine(d_min,d_max):
    for d in range(d_min,d_max+2,2):
        print(find_largest_orbit(discrete_sine_poly(d),int((d+5)/2),int((d+5)/2)))

def create_henon_graphics_discrete_sine(d_min,d_max):
    for d in range(d_min,d_max+2,2):     
        create_henon_graphic(discrete_sine_poly(d),int((d+5)/2),int((d+5)/2),figure_name="outputs/Henon_d_"+str(d),reference_box_size=int((d+5)/2))

def make_your_own_function(values,index_start):
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            print("Outside function range!")
            return None
        return values[x-index_start]
    return p

#poly = make_your_own_function([-10, -2, 2, 3, 2, 0, -2, -3, -2, 2, 10],-5)
#poly = make_your_own_function([-298, -1, 1, -1, 1, 1, 1, 1, 0, -1, -1, -1, -1, 1, -1, 1, 298],-8)

print_largest_orbits_discrete_sine(3,49)

#create_henon_graphics_discrete_sine(3,50) # Runs the original program that was included.