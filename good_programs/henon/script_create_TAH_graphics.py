# Create all the Henon graphics of a certain type (e.g. longest) for the TAH map for each degree d in a range

from henon_tools import *

### Change variables as desired

d_min = 3
d_max = 9

colour_style = "LENGTH"
figure_name_prefix = "TAH_"+colour_style+"_d_"

###

for d in range(d_min,d_max+2,2):     
        create_henon_graphic(discrete_sine_poly(d),int((d+5)/2),
                             figure_name=figure_name_prefix+str(d),
                             colour_style=colour_style)