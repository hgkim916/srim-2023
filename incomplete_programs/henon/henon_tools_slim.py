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
                if len(orbit) == longest_cycle:
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
    plt.savefig(figure_name)
    plt.close()

def make_your_own_function(values,index_start): # define a function that takes specified values at integer intervals starting at index_start.
    def p(x):
        if x-index_start<0 or x-index_start>=len(values):
            #print("Outside function range!")
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

def count_preper(p,radius,x_coefficient):      # returns the total number of preperiodic points of the Henon map of polynomial p in the box given by range
    found_points = []

    for x_tocheck in range(-radius,radius):
        for y_tocheck in range(-radius,radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],radius,x_coefficient=x_coefficient)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
    found_points = [list(tupl) for tupl in {tuple(item) for item in found_points }]   # sanity check, ensures no duplicates     
    return len(found_points)

def shift_poly_in_x(shift,poly):        # shifts the value taken by distance "shift" to the right
                                    # outputs poly(x-shift)
    def new_poly(x):
        return poly(x-shift)
    return new_poly

def next_list_lexicographic_order(input_list,min_entry,max_entry): # returns the next list by alphabetic order
    if min(input_list) >= max_entry:
        return
    elif input_list[0] == max_entry:
        return [min_entry]+next_list_lexicographic_order(input_list[1:],min_entry,max_entry)
    else:
        input_list[0] += 1
        return input_list

def list_of_differences(input_list):
    new_list = [0]*(len(input_list)-1)
    for i in range(len(new_list)):
        new_list[i] = input_list[i+1]-input_list[i]
    return new_list

def find_degree_of_interpolating_polynomial(list_of_values):
    if list_of_values == []: return 0
    if max(list_of_values) == min(list_of_values):
        return 0
    else:
        return 1 + find_degree_of_interpolating_polynomial(list_of_differences(list_of_values))

def print_longest_cycles(f_min,f_max,x_min,x_max,desired_degree):
    f_range = f_max-f_min+1
    x_range = x_max-x_min+1
    a = [f_min]*(x_range)

    biggest = -1

    count = 0
    count2 = 0

    print("range of f:",f_min,f_max)
    print("range of x:",x_min,x_max)
    print("desired degree: <=",desired_degree)
    print("printing max cycle found as we go")
    print("----------")
    print("looking at "+str((f_range)**(x_range))+" values")

    while a != None:
        count2 += 1
        if count2 % 1000000 == 0:
            print(count2)
        poly = make_your_own_function(a,x_min)
        
        if find_degree_of_interpolating_polynomial(a) > desired_degree:
            a = next_list_lexicographic_order(a,f_min,f_max)
            continue

        current = find_longest_cycle_length(poly,max(abs(x_min),abs(x_max)),max(abs(x_min),abs(x_max)))
        if current >= biggest:
            biggest=current
            count += 1
            print(count,a,current,find_degree_of_interpolating_polynomial(a))
        a = next_list_lexicographic_order(a,f_min,f_max)
        
        
    print(count)

def new_family_poly(d):
    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly2(d): # for odd d
    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    list_of_values[-2] = -1
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly3(d): # for even d
    list_of_values = [0]*(d+1)
    list_of_values[0] = 1
    list_of_values[d//2] = 1
    list_of_values[d//2 + 1] = 1
    return make_your_own_function(list_of_values,-(d//2))
#for d in range(3,51,2):
#    roughly_half = (d-1)//2
#    print(find_longest_cycle_length(new_family_poly(d),roughly_half+1,roughly_half+1))


#print_longest_cycles(-1,1,-4,5,desired_degree=9)
#d = 11
#roughly_half = (d-1)//2
#create_henon_graphic(new_family_poly(d),roughly_half+1,roughly_half+1,colour_style="LENGTH")
d = 29
poly = shift_poly_in_x(2,discrete_sine_poly(d))
create_henon_graphic(poly,d,d,colour_style="LONGEST",figure_size=10)