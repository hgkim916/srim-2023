# Create all the Henon graphics of a certain type (e.g. longest) for the specified new family for each degree d in a range

from henon_tools import *

def new_family_poly1(d):
    '''
    Returns the polynomial of degree d in New Family 1.

    Parameters:
        d (int): The degree of the polynomial - must be odd.
    '''
    if d%2 == 0:
        raise ValueError("d must be odd.")

    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly2(d): # for odd d
    '''
    Returns the polynomial of degree d in New Family 2, for odd d.

    Parameters:
        d (int): The degree of the polynomial - must be odd.
    '''
    if d%2 == 0:
        raise ValueError("d must be odd.")
    
    roughly_half = (d-1)//2
    list_of_values = ([1]*(roughly_half)+[0])*2
    list_of_values[-2] = -1
    #print(list_of_values)
    return make_your_own_function(list_of_values,-roughly_half)

def new_family_poly3(d): # for even d
    '''
    Returns the polynomial of degree d in New Family 3, for even d.

    Parameters:
        d (int): The degree of the polynomial - must be even.
    '''
    if d%2 == 1:
        raise ValueError("d must be even.")

    list_of_values = [0]*(d+1)
    list_of_values[0] = 1
    list_of_values[d//2] = 1
    list_of_values[d//2 + 1] = 1
    return make_your_own_function(list_of_values,-(d//2))

def new_family_poly(d,family_index):
    '''
    Returns the polynomial of degree d in the New Family given by family_index.
    
    Parameters:
        d (int): The degree of the polynomial.
        family_index (int): The index of the New Family. Should be 1, 2 or 3.
    '''
    if family_index == 1:
        return new_family_poly1(d)
    elif family_index == 2:
        return new_family_poly2(d)
    elif family_index == 3:
        return new_family_poly3(d)
    else:
        print("family_index must be 1, 2 or 3")
        return None

### Change variables as desired

family = 1          # choices: 1,2,3.  Note family 1 and 2 are defined on odd d, family 3 is defined on even d. 

# Range of d
d_min = 3           # Family 2 requires d to be odd, family 3 requires d to be even.
d_max = 9           # if only one value of d is needed, make d_min=d_max = required value

colour_style = "LENGTH"

figure_name_prefix = "new_family_"+str(family)+"_"+colour_style+"_d_" # Names the output files

###

for d in range(d_min,d_max+2,2):     
        create_henon_graphic(new_family_poly(d,family),
                             radius=(d+1)//2,
                             figure_name=figure_name_prefix+str(d),
                             colour_style=colour_style)