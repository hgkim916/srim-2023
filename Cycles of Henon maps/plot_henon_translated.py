#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

def poly(d):            # returns the polynomial p
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
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
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
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
    elif d%4==2: # given by -1 -1 0 1 1 0 -1 ... (starting at -1/2)
        def p(x):
            x = int(x+1/2)
            # define it periodically 
            if x%3 == 2:
                val = 0
            elif x%6 == 3 or x%6 == 4:
                val = 1
            else:
                val = -1
            x=x-1/2
            # adjust the endpoints and return the appropriate value of p(x)
            if abs(x)<int(d/2)+1:
                return val
            elif abs(x)<int(d/2)+2:
                return val+1
            elif abs(x)<int(d/2)+3:
                return val + (d+2)
            else:  # can only calculate for x in the "nice" range
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None        # This is a flag telling us that we went outside of the nice box, and hence we should stop iterating
    elif d%4==0: # given by 1 1 0 -1 -1 0 1 ... (starting at -1/2)
        def p(x):
            x = int(x+1/2)
            # define it periodically 
            if x%3 == 2:
                val = 0
            elif x%6 == 3 or x%6 == 4:
                val = -1
            else:
                val = 1
            x=x-1/2
            # adjust the endpoints and return the appropriate value of p(x)
            if abs(x)<int(d/2)+1:
                return val
            elif abs(x)<int(d/2)+2:
                return val+1
            elif abs(x)<int(d/2)+3:
                return val + (d+2)
            else:  # can only calculate for x in the "nice" range
                # print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None        # This is a flag telling us that we went outside of the nice box, and hence we should stop iterating
    else:
        print("Warning, this value of d has not yet been implemented!")
        def p(x):       # returns the identity, so that the program doesn't crash
            return x
    return p

def translate_to_dyn_compr(d):
    if d%2==1:
        shift = (d+7)/2
        init_pol = poly(d)
        def p(x):
            if init_pol(x - shift) == None:
                return None                 # carry over the flag
            else:
                return init_pol(x - shift) + 2*shift
        return p 
    else:
        shift = (d+7)/2
        init_pol = poly(d)
        def p(x):
            if init_pol(x - shift) == None:
                return None                 # carry over the flag
            else:
                return init_pol(x - shift) + 2
        return p 
    
def Henon(p,X):     # this is the Henon map of polynomial p: (x,y) -> (y, -x + p(y))
    poly_val = p(X[1])
    if poly_val == None:
        return None         # this is a flag, telling that we went outside of the box
    else:
        return [X[1],-X[0]+poly_val]

def trace_pt(p,X,box_range):        # follows the orbit of a point X under Henon map with polynomial p, 
                                    # stopping either upon repeating a vertex or when we go outside the nice range
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = Henon(p,X)
        if X == None:
            return []               # if we iterate outside the box, we don't want to plot the orbit at all
    orbit.append(X)                 # add the last vertex to complete the orbit
    return orbit

def plot_orbit (orbit,col_index,box_range):         # plot the orbit
    # initialise the colours
    colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
    colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
    col = colours[col_index % len(colours)]         # pick the corresponding colour


    if len(orbit) == 2:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=300/box_range,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit)-1)]
        ys = [orbit[i][1] for i in range(len(orbit)-1)]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices

        for k in range(len(orbit)-1):
            plt.arrow(orbit[k][0],orbit[k][1],orbit[k+1][0]-orbit[k][0],orbit[k+1][1]-orbit[k][1],width=.01,color = col,alpha =0.2)         # plots an arrow between two consecutive iterates

for d in range(1,50):     # save all of these figures
    p = translate_to_dyn_compr(d)     # get the polynomial
    box_range = d+6
    check_range = d+6      # could start cycles in a subset of the box 
    plotted = []        # store all the vertices whose orbits we've already plotted

    col_index = 0
    plt.figure(figsize = (15,15))

    xs = [1,1,box_range,box_range,1]
    ys = [1,box_range,box_range,1,1]
    plt.plot(xs,ys,"--",color = "grey")     # draw the box where p is "nice" to show where the cycles are relative to it
    # plt.grid(which="both")
    # plt.xticks([i for i in range(-box_range-1,box_range+1)])
    # plt.yticks([i for i in range(-box_range-1,box_range+1)])

    for i in range(1,check_range+1):    
        for j in range(1,check_range+1):        # iterate through all the points
            if [i,j] not in plotted:            # only iterate if we haven't plotted already, to reduce computation
                orbit = trace_pt(p,[i,j],box_range)     # get the orbit by tracing the point
                for pt in orbit:
                    plotted.append(pt)          # add each iterate to the list of plotted vertices
                plot_orbit(orbit,col_index,box_range)       # plot the orbit    
                col_index +=1                   # to vary our colours in a nice way

    # plot formatting:

    plt.axis('equal')
    plt.savefig("Cycles_of_Henon_maps/transl_plots/Henon_cyc_d="+str(d)+"")
    plt.close()