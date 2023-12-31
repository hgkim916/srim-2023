#!/usr/bin/env python3


# THE AIM OF THIS:  given any d, it checks which of the following four maps
#                       1. (x,y)->(y,-x+p(y))
#                       2. (x,y)->(y, x+p(y))
#                       3. (x,y)->(y,-x-p(y))
#                       4. (x,y)->(y, x-p(y))
#                   is the best one in terms of the maximum box which must be periodic.


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
            if abs(x)<=d:
                return val
            elif x == d+1:
                return val+1
            elif x == -(d+1):
                return val-1
            elif x == d+2:
                return val + (d+2)
            elif x == -(d+2):
                return val - (d+2)
            else:  # can only calculate for x in the "nice" range
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return 0
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
            if abs(x)<=d:
                return val
            elif x == d+1:
                return val+1
            elif x == -(d+1):
                return val-1
            elif x == d+2:
                return val + (d+2)
            elif x == -(d+2):
                return val - (d+2)
            else: # can only calculate for x in the "nice" range
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return 0
    else:
        print("Warning, this value of d has not yet been implemented!")
        def p(x):       # returns the identity, so that the program doesn't crash
            return x
        
    return p
            
def Henon(p,X,i):     # this is the Henon map of polynomial p: (x,y) -> (y, -x + p(y))
    if i == 1:
        return [X[1],-X[0]+p(X[1])]
    elif i == 2:
        return [X[1], X[0]+p(X[1])]
    elif i == 3:
        return [X[1],-X[0]-p(X[1])]
    elif i == 4:
        return [X[1], X[0]-p(X[1])]
    else:
        print("ERROR: Henon map with i=",i," is undefined!")
        return X

def trace_pt(p,X,box_range,henon_index):        # follows the orbit of a point X under Henon map with polynomial p, 
                                    # stopping either upon repeating a vertex or when we go outside the nice range
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = Henon(p,X,henon_index)
        if abs(X[0])>box_range or abs(X[1])>box_range:
            return []               # if we iterate outside the d+2 box, we don't want to plot the orbit at all
    orbit.append(X)                 # add the last vertex to complete the orbit
    return orbit

def find_size(d,visited,henon_index):       # prints the size of the max l_inf box allowed
    for size in range(0,d+3):
            for i in range(-size,size+1):
                if [i,size] not in visited or [i,-size] not in visited or [size,i] not in visited or [-size,i] not in visited:
                    print("Maximum size for the",henon_index,"th Henon map for d=",d," is d -",d-(size-1))
                    return 
                

def output_bound():         # OUTPUT BOUND
    for d in range(11,40,12):     
        p = poly(d)     # get the polynomial
        box_range = d+2
        check_range = d+2       # could start cycles in a subset of the box 

        for hen_index in range(1,5):
            visited = []        # store all the vertices whose orbits we've already plotted
            
            for i in range(0,check_range+1):    
                for j in range(0,check_range+1):        # iterate through all the points
                    if [i,j] not in visited:            # only iterate if we haven't plotted already, to reduce computation
                        orbit = trace_pt(p,[i,j],box_range,hen_index)     # get the orbit by tracing the point
                        for pt in orbit:
                            visited.append(pt)          # add each iterate to the list of plotted vertices
            find_size(d,visited,hen_index)


output_bound()