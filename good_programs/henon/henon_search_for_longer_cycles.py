from henon_tools_slim import *

def next_list_colex_order(input_list,min_entry,max_entry): 
    '''
    Returns the next list by colexicographic order, given input_list, the minimum value an entry should take,
    and the maximum value an entry should take.

    This is like, the next list alphabetically, in a dictionary, if everyone read right to left.

    For a demo see demo_lexicographic_lists.py
    '''
    if min(input_list) >= max_entry:
        return
    elif input_list[0] == max_entry:
        return [min_entry]+next_list_colex_order(input_list[1:],min_entry,max_entry)
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
            a = next_list_colex_order(a,f_min,f_max)
            continue

        current = find_longest_cycle_length(poly,max(abs(x_min),abs(x_max)),max(abs(x_min),abs(x_max)))
        if current >= biggest:
            biggest=current
            count += 1
            print(count,a,current,find_degree_of_interpolating_polynomial(a))
        a = next_list_colex_order(a,f_min,f_max)
        
    print(count)