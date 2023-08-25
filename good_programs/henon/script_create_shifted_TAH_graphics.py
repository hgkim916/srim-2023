# Create all the Henon graphics of a certain type (e.g. longest) for a specified shift of the TAH polys for each degree d in a range

from henon_tools import *

def shift_poly_in_x(shift,poly):       
    '''
    Shifts the given function to the right, returning poly(x-shift).

    Parameters:
        shift (int): The distance to shift the function by.
        poly (function): The function to shift.
    '''
    def new_poly(x):
        return poly(x-shift)
    return new_poly

### Change variables as desired

# Range of d
d_min = 3
d_max = 9
step = 6

shift = 2 # Choose which shift to plot

colour_style = "LENGTH"
figure_name_prefix = "shifted_by_"+str(shift)+"_TAH_"+colour_style+"_d_"

###

for d in range(d_min,d_max+2,step):     
        create_henon_graphic(shift_poly_in_x(shift,discrete_sine_poly(d)),int((d+5)/2)+shift,
                             figure_name=figure_name_prefix+str(d),
                             colour_style=colour_style)