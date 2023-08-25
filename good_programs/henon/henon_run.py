from henon_tools_slim import *


'''
EXAMPLES OF WHAT TO RUN:

1. Create a single Henon graphic
    a. TAH
    b. Shifted TAH
    c. New Family i
2. Create all the Henon graphic for each degree d in a range


------------------------------------------------------------------------------------------------------------------------

1. Create a single Henon graphic
    a. TAH:
'''
# Change as desired:
d = 19
chosen_colour_style = "LENGTH"
filename = "TAH_d"+str(d)+"_graphic"

# Remove '# ' from the line below to run the function and create the Henon graphic
# create_henon_graphic(discrete_sine_poly(d),escape_radius=(d+5)//2,check_radius=(d+5)//2,figure_name=filename,colour_style=chosen_colour_style,figure_size=10)


'''
    b. Shifted TAH:
'''
# Change as desired:
d = 19
shift = 1
chosen_colour_style = "LENGTH"
filename = "Shifted_TAH_d"+str(d)+"_graphic"

# Remove '# ' from the line below to run the function and create the Henon graphic
create_henon_graphic(shift_poly_in_x(shift,discrete_sine_poly(d)),escape_radius=(d+5)//2,check_radius=(d+5)//2,colour_style=chosen_colour_style,figure_size=10)

'''
    c. New Family i:
'''
# Change the degree d, family_index and colour style as desired
d = 19
family_index = 1
chosen_colour_style = "LENGTH"
filename = "New_Family_"+str(family_index)+"_d"+str(d)+"_graphic"

# Remove '# ' from the line below to run the function and create the Henon graphic
# create_henon_graphic(new_family_poly(d,family_index),escape_radius=(d+1)//2,check_radius=(d+1)//2,colour_style=chosen_colour_style,figure_size=10)
