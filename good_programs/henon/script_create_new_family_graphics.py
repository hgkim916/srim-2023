# Create all the Henon graphics of a certain type (e.g. longest) for the specified new family for each degree d in a range

from henon_tools_slim import *

# Change as desired:

family = 1          # choices: 1,2,3.  

d_min = 3           # Family 2 requires d to be odd, family 3 requires d to be even.
d_max = 9           # if only one value of d is needed, make d_min=d_max = required value
colour_style = "LENGTH"
figure_name_prefix = "new_family_"+str(family)+"_"+colour_style+"_d_" # Saves as e.g. "TAH_LENGTH_d_19.png" 

for d in range(d_min,d_max+2,2):     
        create_henon_graphic(new_family_poly(d,family),escape_radius=(d+1)//2,check_radius=(d+1)//2,
                             figure_name=figure_name_prefix+str(d),
                             colour_style=colour_style)