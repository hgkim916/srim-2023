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

def plot_orbit (orbit,colour_parameter,box_range,colour_style="DEFAULT"):         # plot the orbit
    
    if colour_style == "DEFAULT":
        # initialise the colours
        colours=["navy","mediumblue","slateblue","blueviolet","indigo","mediumorchid","thistle","plum","magenta","deeppink","crimson","lightpink","salmon","red","brown","maroon","saddlebrown","peru","sandybrown","lightsalmon","darkorange","goldenrod","gold","khaki","y","olive","olivedrab","yellowgreen","chartreuse", "limegreen", "g", "seagreen","mediumaquamarine","lightseagreen","teal","c","aqua","deepskyblue","lightskyblue","steelblue"]
        colours = [colours[3*i%len(colours)] for i in range(0,len(colours))] # shuffle so that they're not too close
        col = colours[colour_parameter % len(colours)]         # pick the corresponding colour
    if colour_style == "PARAMETER":
        col = cm.cool(colour_parameter)

    if len(orbit) == 1:     # this must be a fixed point
        plt.plot(orbit[0][0],orbit[0][1],marker=r'$\circlearrowleft$',ms=300/box_range,color = col)
    else:
        xs = [orbit[i][0] for i in range(len(orbit))]
        ys = [orbit[i][1] for i in range(len(orbit))]
        plt.scatter(xs,ys,color = col)              # plots the individual vertices
        orbit.append(orbit[0])      # adds the first vertex to the end, to create a loop
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
                if len(orbit)> largest_orbit_size:
                    #print(orbit)
                    largest_orbit_size = len(orbit)
    return largest_orbit_size

def create_henon_graphic(p,escape_radius,check_radius,figure_name="output",reference_box_size=0,colour_style="DEFAULT",x_coefficient=-1):

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
        longest_cycle = find_longest_cycle_length(p,escape_radius,check_radius,x_coefficient=x_coefficient)
    if colour_style == "LONGEST":
        exceptional_cycle_already_plotted = False

    count = 0
    for i in range(-check_radius,check_radius+1):    
        for j in range(-check_radius,check_radius+1):        # iterate through all the points
            if [i,j] in found_points: continue            # only iterate if we haven't plotted already, to reduce computation
            
            orbit = trace_pt(p,[i,j],escape_radius,x_coefficient=x_coefficient)     # get the orbit by tracing the point
            if len(orbit) == 0: continue
            print(orbit)
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

def print_length_of_longest_cycles_discrete_sine(d_min,d_max,negative = False,x_coefficient=-1):
    for d in range(d_min,d_max+2,2):
        print(find_longest_cycle_length(discrete_sine_poly(d,negative=negative),int((d+5)/2),int((d+5)/2),x_coefficient=x_coefficient))

def create_henon_graphics_discrete_sine(d_min,d_max,figure_name="henon_d_",colour_style = "DEFAULT",negative = False,x_coefficient=-1):
    for d in range(d_min,d_max+2,2):     
        create_henon_graphic(discrete_sine_poly(d,negative=negative),int((d+5)/2),int((d+5)/2),
                             figure_name="outputs/"+figure_name+str(d),
                             reference_box_size=int((d+5)/2),
                             colour_style=colour_style,
                             x_coefficient=x_coefficient)

def make_your_own_function(values,index_start): # define a function that takes specified values at integer intervals starting at index_start.
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            print("Outside function range!")
            return None
        return values[x-index_start]
    return p

def count_cycle_lengths(p,escape_radius,check_radius,x_coefficient=-1): # returns a dict with all the counts of cycle lengths
    found_points = []
    lengths = {}
    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius,x_coefficient=x_coefficient)
                orbit_length = len(orbit)
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

def count_preper(p,range,x_coefficient):      # returns the total number of preperiodic points of the Henon map of polynomial p in the box given by range
    found_points = []

    for x_tocheck in range(-range,range):
        for y_tocheck in range(-range,range):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],range,x_coefficient=x_coefficient)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
                
    return len(found_points)

def shift_poly_in_x(shift,poly):        # shifts the value taken by distance "shift" to the right
    def new_poly(x):
        return poly(x-shift)
    return new_poly

def check_cycles_left_right_shifts(dmin,dmax,step):
    for d in range(dmin,dmax,step):
        f = open("longest_cycle_shifts.txt",'a')
        init_poly = discrete_sine_poly(d)
        shift_min = -4
        shift_max = 4
        max_cycle_lengths = []
        for shift in range(shift_min,shift_max):
            shift_poly = shift_poly_in_x(shift,init_poly)
            max_cycle_lengths.append(find_longest_cycle_length(shift_poly,int((d+5)/2)+abs(shift),int((d+5)/2)+abs(shift)))
            
        max_length = max(max_cycle_lengths)
        print("d=",d," - longest cycle:",max_length,"     achieved by shifting by: ",end="",file=f)
        for shift in range(shift_min,shift_max):
            if max_cycle_lengths[shift-shift_min] == max_length:
                print(shift,end=" ",file=f)
        print("",file=f)
        f.close()

def check_no_of_preper_shifted_polys(dmin,dmax,step):
    for d in range(dmin,dmax,step):
        f = open("no_of_preper_of_shifts.txt",'a')
        print("d=",d,":",file=f)
        f.close()
        init_poly = discrete_sine_poly(d)
        shift_min = -2
        shift_max = 2
        for shift in range(shift_min,shift_max+1):
            shift_poly = shift_poly_in_x(shift,init_poly)
            f = open("no_of_preper_of_shifts.txt",'a')
            print("shift =",shift,"gives",count_preper(shift_poly,int((d+1)/2)+abs(shift),-1),"preperiodic pts",file=f)
            f.close()
        f = open("no_of_preper_of_shifts.txt",'a')
        print(" ",file=f)
        f.close()

def check_no_of_preper_AND_max_cycle_shifted_polys(dmin,dmax,step):
    for d in range(dmin,dmax,step):
        f = open("no_of_preper_and_max_cycle_of_shifts.txt",'a')
        print("d=",d,":",file=f)
        f.close()
        init_poly = discrete_sine_poly(d)
        shift_min = -2
        shift_max = 2
        for shift in range(shift_min,shift_max+1):
            shift_poly = shift_poly_in_x(shift,init_poly)
            f = open("no_of_preper_and_max_cycle_of_shifts.txt",'a')
            print("shift =",shift,"gives",count_preper(shift_poly,int((d+1)/2)+abs(shift),-1),"preperiodic pts and a max cycle of length",find_longest_cycle_length(shift_poly,int((d+5)/2)+abs(shift),int((d+5)/2)+abs(shift)),file=f)
            f.close()
        f = open("no_of_preper_and_max_cycle_of_shifts.txt",'a')
        print(" ",file=f)
        f.close()


def create_all_henon_graphics_discrete_sine():
    create_henon_graphics_discrete_sine(3,49,figure_name="default/henon_d_")
    create_henon_graphics_discrete_sine(3,49,figure_name="by_length/henon_bylength_d_",colour_style="LENGTH")
    create_henon_graphics_discrete_sine(3,49,figure_name="longest/henon_longest_d_",colour_style="LONGEST")
    create_henon_graphics_discrete_sine(3,49,figure_name="default/henon_negative_d_",negative=True)
    create_henon_graphics_discrete_sine(3,49,figure_name="by_length/henon_negative_bylength_d_",colour_style="LENGTH",negative=True)
    create_henon_graphics_discrete_sine(3,49,figure_name="longest/henon_negative_longest_d_",colour_style="LONGEST",negative=True)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/positive_discrete_sine/default/henon_d_",x_coefficient=1)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/positive_discrete_sine/by_length/henon_bylength_d_",colour_style="LENGTH",x_coefficient=1)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/positive_discrete_sine/longest/henon_longest_d_",colour_style="LONGEST",x_coefficient=1)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/negative/default/henon_d_",negative=True,x_coefficient=1)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/negative/by_length/henon_bylength_d_",colour_style="LENGTH",negative=True,x_coefficient=1)
    create_henon_graphics_discrete_sine(3,29,figure_name="x_coefficient_1/negative/longest/henon_longest_d_",colour_style="LONGEST",negative=True,x_coefficient=1)


check_cycles_left_right_shifts(3,200,2)
                

#d = 47
#poly = discrete_sine_poly(d) # Check values of the discrete sine
#print(find_longest_cycle_length(poly,int((d+5)/2),int((d+5)/2)))
#create_henon_graphic(poly,int((d+5)/2),int((d+5)/2),figure_name="output",reference_box_size=0,colour_style="LONGEST")

#for i in range(-8,9):
#    print(i,poly(i))
