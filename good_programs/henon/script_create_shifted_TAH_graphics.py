# This will create all the Henon graphics of a certain type (e.g. longest) and a specified shift for each degree d in a range

from henon_tools_slim import *

# Change as desired:

d_min = 3
d_max = 9           # if only one value of d is needed, make d_min=d_max = required value
shift = 2
step = 6            # by default, only looks at d of the same value mod 6, but can change to any other EVEN step
colour_style = "LENGTH"
figure_name_prefix = "shifted_by_"+str(shift)+"_TAH_"+colour_style+"_d_" # Saves as e.g. "TAH_LENGTH_d_19.png" 

for d in range(d_min,d_max+2,step):     
        create_henon_graphic(shift_poly_in_x(shift,discrete_sine_poly(d)),int((d+5)/2)+shift,int((d+5)/2)+shift,
                             figure_name=figure_name_prefix+str(d),
                             colour_style=colour_style)