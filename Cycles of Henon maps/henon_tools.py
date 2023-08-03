#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

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
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
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
                print("Error, cannot calculate for d=",d,"and x=",x,"!")
                return None
            
            if negative:
                val = -val
            return val
        return p
    
    print("Warning, this value of d has not yet been implemented!")
    return None
            
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

def plot_orbit (orbit,colour_parameter,box_range,colour_style = "DEFAULT"):         # plot the orbit
    
    if colour_style == "DEFAULT":
        # initialise the colours
        colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
        colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
        col = colours[colour_parameter % len(colours)]         # pick the corresponding colour
    if colour_style == "PARAMETER":
        col = cm.cool(colour_parameter)

    if len(orbit) == 2:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=300/box_range,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit)-1)]
        ys = [orbit[i][1] for i in range(len(orbit)-1)]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices

        for k in range(len(orbit)-1):
            plt.arrow(orbit[k][0],orbit[k][1],orbit[k+1][0]-orbit[k][0],orbit[k+1][1]-orbit[k][1],width=.01,color = col,alpha =0.2)         # plots an arrow between two consecutive iterates

def find_longest_cycle_length(p,escape_radius,check_radius): # find the length of the longest cycle in an orbit p.
                                                             # can easily be modified to give you the actual cycle too.
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

def create_henon_graphic(p,escape_radius,check_radius,figure_name="output",reference_box_size=0,colour_style="DEFAULT"):

    found_points = []        # store all the vertices whose orbits we've already plotted
    plt.figure(figsize = (15,15))

    if reference_box_size>0: # Draw the reference box
        xs = [-reference_box_size,-reference_box_size,reference_box_size,reference_box_size,-reference_box_size]     
        ys = [-reference_box_size,reference_box_size,reference_box_size,-reference_box_size,-reference_box_size]
        plt.plot(xs,ys,"--",color = "grey")
    # plt.grid(which="both")
    # plt.xticks([i for i in range(-box_range-1,box_range+1)])
    # plt.yticks([i for i in range(-box_range-1,box_range+1)])

    if colour_style == "LENGTH" or "LONGEST":
        longest_cycle = find_longest_cycle_length(p,escape_radius,check_radius)
    if colour_style == "LONGEST":
        exceptional_cycle_already_plotted = False

    count = 0
    for i in range(-check_radius,check_radius+1):    
        for j in range(-check_radius,check_radius+1):        # iterate through all the points
            if [i,j] in found_points: continue            # only iterate if we haven't plotted already, to reduce computation
            
            orbit = trace_pt(p,[i,j],escape_radius)     # get the orbit by tracing the point
            if len(orbit) == 0: continue
            found_points.extend(orbit)          # add each iterate to the list of plotted vertices
            if colour_style == "DEFAULT":
                plot_orbit(orbit,count,escape_radius)       # plot the orbit    
            elif colour_style == "LENGTH":
                colour_param = round((len(orbit)-1)/longest_cycle*255)
                plot_orbit(orbit,colour_param,escape_radius,colour_style="PARAMETER")
            elif colour_style == "LONGEST":
                if len(orbit)-1 == longest_cycle:
                    if exceptional_cycle_already_plotted:
                        colour_param = 127
                    else:
                        colour_param = 255
                        exceptional_cycle_already_plotted = True
                else:
                    colour_param = 0
                plot_orbit(orbit,colour_param,escape_radius,colour_style="PARAMETER")
            count += 1                   # for count-based colour styles

    # plot formatting:

    plt.axis('equal')
    plt.savefig(figure_name)
    plt.close()

def print_length_of_longest_cycles_discrete_sine(d_min,d_max,negative = False):
    for d in range(d_min,d_max+2,2):
        print(find_longest_cycle_length(discrete_sine_poly(d,negative=negative),int((d+5)/2),int((d+5)/2)))

def create_henon_graphics_discrete_sine(d_min,d_max,figure_name="henon_d_",colour_style = "DEFAULT",negative = False):
    for d in range(d_min,d_max+2,2):     
        create_henon_graphic(discrete_sine_poly(d,negative=negative),int((d+5)/2),int((d+5)/2),figure_name="outputs/"+figure_name+str(d),reference_box_size=int((d+5)/2),colour_style=colour_style)

def make_your_own_function(values,index_start): # define a function that takes specified values at integer intervals starting at index_start.
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            print("Outside function range!")
            return None
        return values[x-index_start]
    return p

def count_cycle_lengths(p,escape_radius,check_radius): # returns a dict with all the counts of cycle lengths
    found_points = []
    lengths = {}
    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius)
                orbit_length = len(orbit)-1
                if orbit_length == -1:
                    continue
                #print(orbit)
                if orbit_length not in lengths:
                    lengths[orbit_length] = 1
                else:
                    lengths[orbit_length] += 1
                
                found_points.extend(orbit)
    lengths = {a:b for a,b in sorted(lengths.items())}
    return lengths



#d = 47
#poly = discrete_sine_poly(d) # Check values of the discrete sine
#print(find_longest_cycle_length(poly,int((d+5)/2),int((d+5)/2)))
#create_henon_graphic(poly,int((d+5)/2),int((d+5)/2),figure_name="output",reference_box_size=0,colour_style="LONGEST")

#for i in range(-8,9):
#    print(i,poly(i))

#create_henon_graphics_discrete_sine(3,49,figure_name="default/henon_d_")
#create_henon_graphics_discrete_sine(3,49,figure_name="by_length/henon_bylength_d_",colour_style="LENGTH")
#create_henon_graphics_discrete_sine(3,49,figure_name="longest/henon_longest_d_",colour_style="LONGEST")
#create_henon_graphics_discrete_sine(3,49,figure_name="default/henon_negative_d_",negative=True)
#create_henon_graphics_discrete_sine(3,49,figure_name="by_length/henon_negative_bylength_d_",colour_style="LENGTH",negative=True)
#create_henon_graphics_discrete_sine(3,49,figure_name="longest/henon_negative_longest_d_",colour_style="LONGEST",negative=True)

print_length_of_longest_cycles_discrete_sine(3,49,negative=True)