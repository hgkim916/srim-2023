#!/usr/bin/env python3


# THE AIM OF THIS:  given any d, it outputs the maximum size of a box within which all points are stuck in a cycle


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

def trace_pt(p,X):        # follows the orbit of a point X under Henon map with polynomial p, 
                                    # stopping either upon repeating a vertex or when we go outside the nice range
    orbit = []
    while X not in orbit:
        orbit.append(X)
        X = Henon(p,X)
        if X == None:
            return []               # if we iterate outside the box, we don't want to plot the orbit at all
    orbit.append(X)                 # add the last vertex to complete the orbit
    return orbit

def expected_size(d):
    if d%6==1:
        return d+2
    elif d%6==3:
        return d-2
    elif d%6 == 5:
        return d-4
    elif d%2 == 0:
        return 0
    else:
        print("Uhhh I shouldn't have d=",d)
        return 0

def find_size(d,visited):       # prints the size of the max l_inf box allowed
    # find lower corner:
    for i in range(1,int((d+6)/2)+2):
        complete = True
        j = (d+7)-i
        for k in range(i,j+1):
            if [i,k] not in visited or [k,i] not in visited or [j,k] not in visited or [k,j] not in visited:
                complete = False
                break
        if complete:
            upper = j
            lower = i
            return upper - lower + 1
    if i<j:
        return 0

def output_bound():         # OUTPUT BOUND
    for d in range(2,50,2):     
        p = translate_to_dyn_compr(d)     # get the polynomial
        check_range = d+6      # could start cycles in a subset of the box 
        visited = []        # store all the vertices whose orbits we've already plotted

        for i in range(1,check_range+1):    
            for j in range(1,check_range+1):        # iterate through all the points
                if [i,j] not in visited:            # only iterate if we haven't plotted already, to reduce computation
                    orbit = trace_pt(p,[i,j])     # get the orbit by tracing the point
                    for pt in orbit:
                        visited.append(pt)          # add each iterate to the list of plotted vertices
        max_size = find_size(d,visited)
        # print(max_size)
        print("Maximum size for d=",d,"is d + 6 -",d + 6- max_size)


def check_bound():          # CHECK BOUND
    f = open("Cycles_of_Henon_maps/check_henon_transl.txt","a") 
    f.close()   
    f = open("Cycles_of_Henon_maps/check_henon_transl.txt","w") 
    f.close()   
    for d in range(1,300):
        f = open("Cycles_of_Henon_maps/check_henon_transl.txt","a")    
        p = translate_to_dyn_compr(d)     # get the polynomial
        check_range = d+6    # could start cycles in a subset of the box 
        visited = []        # store all the vertices whose orbits we've already plotted
        for i in range(1,check_range+1):    
            for j in range(1,check_range+1):        # iterate through all the points
                if [i,j] not in visited:            # only iterate if we haven't plotted already, to reduce computation
                    orbit = trace_pt(p,[i,j])     # get the orbit by tracing the point
                    for pt in orbit:
                        visited.append(pt)          # add each iterate to the list of plotted vertices
        max_size = find_size(d,visited)
        expected = expected_size(d)
        if max_size!= expected:
            print("The box for d=",d,"has size d+6 -",(d+6)-max_size,", different from the expected size of d+6 -",d+6-expected,file=f)
        else:
            print("d=",d,"- good!",file=f)
        f.close()
    

# output_bound()
check_bound()