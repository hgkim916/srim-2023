# This will create all the Henon graphics of a certain type (e.g. longest) for each degree d in a range

from henon_tools_slim import *

# Change as desired:

d_min = 3
d_max = 71
figure_name = "TAH_d_" # Saves as e.g. "TAH_d_19.png" 
colour_style = "DEFAULT"

for d in range(d_min,d_max+2,2):     
        create_henon_graphic(discrete_sine_poly(d),int((d+5)/2),int((d+5)/2),
                             figure_name=figure_name+str(d),
                             reference_box_size=int((d+5)/2),
                             colour_style=colour_style)