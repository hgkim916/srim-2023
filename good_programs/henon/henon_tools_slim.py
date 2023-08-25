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
        x_coeff (int): The delta in the general Henon map. Default is -1.
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

def find_longest_cycle_length(p,escape_radius,check_radius,x_coeff=-1):
    '''
    Returns the length of the longest periodic cycle of the Henon map (x,y) -> (y,x_coeff*x + p(y)) within the box of radius check_radius,.
    Parameters:
        p (function): The polynomial associated with the Henon map
        escape_radius (int): If the point iterates outside the box of radius escape_radius, then assume it escapes.
        check_radius (int): The box to check for periodic cycles.
        x_coeff (int): The delta in the general Henon map. Default is -1.
    '''

    found_points = []
    largest_orbit_size = 0

    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius,x_coeff=x_coeff)
                if len(orbit) == 0:
                    continue
                #print(orbit)
                found_points.extend(orbit)
                if len(orbit) > largest_orbit_size:
                    #print(orbit)
                    largest_orbit_size = len(orbit)
    return largest_orbit_size

def create_henon_graphic(p,escape_radius,check_radius
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
        escape_radius (int): If the point iterates outside the box of radius escape_radius, then assume it escapes.
        check_radius (int): The box to check for periodic cycles.
        figure_name (str): (Default: "output") The name of the file to save the figure to (extension is assumed to be .png)
        figure_title (str): (Default: no title) The title of the figure (shown on the plot)
        figure_size (int): (Default: 10) The size of the figure
        reference_box_size (int): (Default: 0, i.e. no reference box) The size of the reference box to draw
        colour_style (str): (default: "DEFAULT") Options to choose the colour of the plot.
            "DEFAULT": A collection of arbitrary colours.
            "LENGTH": The colour of each cycle is determined by its length, with the longest cycle being magenta and the shortest being blue.
            "LONGEST": The colour of each cycle is determined by whether it is the longest cycle or not.
        x_coeff (int): The delta in the general Henon map. Default is -1.
    '''

    found_points = []        # store all the vertices whose orbits we've already plotted
    fig,ax = plt.subplots(figsize = (figure_size,figure_size))

    if reference_box_size>0: # Draw the reference box
        xs = [-reference_box_size,-reference_box_size,reference_box_size,reference_box_size,-reference_box_size]     
        ys = [-reference_box_size,reference_box_size,reference_box_size,-reference_box_size,-reference_box_size]
        plt.plot(xs,ys,"--",color = "grey")

    if colour_style == "LENGTH" or "LONGEST":
        longest_cycle = find_longest_cycle_length(p,escape_radius,check_radius,x_coeff=x_coeff)
    if colour_style == "LONGEST":
        exceptional_cycle_already_plotted = False

    count = 0
    for i in range(-check_radius,check_radius+1):    
        for j in range(-check_radius,check_radius+1):        # iterate through all the points
            if [i,j] in found_points: continue            # only iterate if we haven't plotted already, to reduce computation
            
            orbit = trace_pt(p,[i,j],escape_radius,x_coeff=x_coeff)     # get the orbit by tracing the point
            if len(orbit) == 0: continue
            #print(orbit)
            found_points.extend(orbit)          # add each iterate to the list of plotted vertices
            if colour_style == "DEFAULT":
                plot_orbit(orbit,count,figure_scale=check_radius/figure_size)       # plot the orbit    
            elif colour_style == "LENGTH":
                colour_param = round((len(orbit)-1)/longest_cycle*255)
                plot_orbit(orbit,colour_param,figure_scale=check_radius/figure_size,colour_style="PARAMETER")
            elif colour_style == "LONGEST":
                if len(orbit) == longest_cycle:
                    if exceptional_cycle_already_plotted:
                        colour_param = 127
                    else:
                        colour_param = 255
                        exceptional_cycle_already_plotted = True
                else:
                    colour_param = 0
                plot_orbit(orbit,colour_param,figure_scale=check_radius/figure_size,colour_style="PARAMETER")
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
    Define a function that takes specified values at integers starting at index_start.
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

def count_cycle_lengths(p,escape_radius,check_radius,x_coeff=-1):
    '''
    Returns a dictionary with all the counts of cycle lengths.
    Parameters: 
        p (function): The polynomial associated with the Henon map
        escape_radius (int): If the point iterates outside the box of radius escape_radius, then assume it escapes.
        check_radius (int): The box to check for periodic cycles.
        x_coeff (int): The delta in the general Henon map. Default is -1.
    '''
    found_points = []
    lengths = {}
    for x_tocheck in range(-check_radius,check_radius):
        for y_tocheck in range(-check_radius,check_radius):
            if [x_tocheck,y_tocheck] in found_points:
                continue
            else:
                orbit = trace_pt(p,[x_tocheck,y_tocheck],escape_radius,x_coeff=x_coeff)
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
    Returns the total number of integer preperiodic points of the Henon map of polynomial p in the box given by radius.
    Parameters:
        p (function): The polynomial associated with the Henon map
        radius (int): The box to check for preperiodic points.
        x_coeff (int): The delta in the general Henon map. Default is -1.
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

def shift_poly_in_x(shift,poly):       
    '''
    Shifts the value taken by distance "shift" to the right, returning poly(x-shift).
    Parameters:
        shift (int): The distance to shift the function by.
        poly (function): The function to shift.
    '''
    def new_poly(x):
        return poly(x-shift)
    return new_poly

def next_list_lexicographic_order(input_list,min_entry,max_entry): 
    '''
    Returns the next list by alphabetic order. VIVIAN PLEASE CLARIFY WHAT THIS IS, I'M A BIT CONFUSED.
    '''
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
    '''
    Prints the longest cycles of the Henon map (x,y) -> (y,x_coeff*x + p(y)) for all polynomials p with degree <= desired_degree.

    Searches amongs polynomials p taking integer values between f_min and f_max on the integers between x_min and x_max.

    Parameters:
        f_min (int): The minimum value that the polynomial can take on the given range.
        f_max (int): The maximum value that the polynomial can take on the given range.
        x_min (int): The first integer to interpolate the function at.
        x_max (int): The last integer to interpolate the function at.
    '''
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
    '''
    Returns the polynomial of degree d in New Family 1

    Parameters:
        d (int): The degree of the polynomial.
    '''
    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly2(d): # for odd d
    '''
    Returns the polynomial of degree d in New Family 2, for odd d.

    Parameters:
        d (int): The degree of the polynomial. Should be odd.
    '''
    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    list_of_values[-2] = -1
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly3(d): # for even d
    '''
    Returns the polynomial of degree d in New Family 3, for even d.

    Parameters:
        d (int): The degree of the polynomial. Should be even.
    '''
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
d = 3
poly = shift_poly_in_x(0,discrete_sine_poly(d))

for i in range((-d-7)//2,(d+7)//2+1):
    print(i,poly(i))

#create_henon_graphic(poly,d,d,colour_style="LONGEST",figure_size=10)