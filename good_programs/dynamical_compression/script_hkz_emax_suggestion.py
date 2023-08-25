from find_hkz_bounds import *

def find_emax_upperbound_suggestion(d, print_progress = True, print_result = True):
    '''
    Finds the "largest" e until HKZ reduction tells us that we can't have
    dynamical compression from [d+e+1] -> [d+e+1]. (i.e. until 4 consecutive failures)
    Note if the 3rd successive minimum has max-norm greater than ceil((d+e)/2), then we can't have dynamical compression.

    Returns the suggested upper bound for e.
    
    Parameters:
        d (int): The degree we want to find a suggestion for. Beyond a certain value (around 55-60), the program gets extremely slow.
        print_progress (bool): (Default: True) Prints the current progress of the function. 
                               Also useful to check for any meaningful floating-point error. 
        print_result (bool): (Default: True)
    '''
    max_e = 10
    e = 1
    while e <= max_e:
        A = lambda_d_e_basis(d,e)
        lower_bound_3 = hkz_succ_min_bounds_maxnorm(A,[3])[0][3]
        if print_progress:
            print(d,e,lower_bound_3,math.floor((d+e+1)/2),lower_bound_3<math.floor((d+e+1)/2))
        if lower_bound_3<math.floor((d+e+1)/2):
            max_e = e + 4 # Require 4 Falses in a row to give up
        e += 1
    if print_result:
        print("----------------",d,max_e - 4)
    return max_e - 4

find_emax_upperbound_suggestion(60,print_progress = True)
