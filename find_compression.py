from fpylll import *
from fpylll.fplll.gso import MatGSO
import math
import polynomial_tools

def lambda_d_e_basis(d,e):
    '''Returns the basis (u_0, u_1, ..., u_(d+1)) as an IntegerMatrix object, with rows as basis vectors.
    '''
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def print_dynamical_compression(degree,compression,number_of_solutions,
                               show_only_dyn_comp = True,
                               show_vector = True,
                               show_polynomial = False,
                               show_translated_poly = False,
                               print_all_failures = True):
    '''
    Prints all possible dynamically compressing IVPs of degree d which send [d+compression] -> [d+compression],
        not repeating identical polynomials up to constants.
    Returns whether a dynamically compressing polynomial was found.
    
    Parameters:
        degree (int): See above
        compression (int): See above
        number_of_solutions (int): The total number of solutions to display. 
                                   This must be finite, due to fpylll limitations.
        show_only_dyn_comp (bool): (Default: True) Whether to only display genuinely dynamically compressing polynomials. 
                                                  If false (default), e.g. a polynomial sending [10] -> [15] would display.
        show_vector (bool): (Default: True) In the output, also print the associated vector in the lattice LambdaDE.
        show_polynomial (bool): (Default: False) In the output, also print the actual polynomial.
        show_translated_poly (bool): (Default: False) In the output, also print the associated dynamically compressing polynomial.
        print_all_failures (bool): (Default: True) If false, hides when the same polynomial up to a constant is found.
                                                   Can cut down on the amount of outsput significantly.
    '''
    # First check for validity of input
    if compression < 1:
        raise ValueError("Invalid compression!")
    
    print("----------PossibleDynCompressing")
    print("Looking for IVPs of degree",degree,"which compress ["+str(degree+compression)+"]")

    # This part is basically taken from fpylll's Github page. We're not exactly sure how it works.
    A = lambda_d_e_basis(degree,compression-1)
    M = MatGSO(A)
    _ = M.update_gso() # What does this do?
    enum = Enumeration(M,strategy = EvaluatorStrategy.BEST_N_SOLUTIONS,nr_solutions = number_of_solutions)
    euclidean_bound = (degree+compression)*(math.floor((degree+compression)/2)**2) # If there exists any compressing polynomials, there is one within this box.
    e1 = enum.enumerate(0, degree+1, euclidean_bound, 0)
    count = 0
    found_vectors = []
    found_dyn_comp = False
    print("Found",len(e1),"possible solutions")
    for i,j in e1:
        count += 1
        # Check for unwanted cases
        if max(j[2:]) == 0 and min(j[2:]) == 0:
            print(count,"(linear solution)")
            continue
        if j[-1] == 0:
            print(count,"(lower degree solution)")
            continue
        #print((i),tuple(map(int,j)))
        
        # Check if we found the same polynomial up to a constant
        non_constant_part = j[1:]
        if non_constant_part in found_vectors:
            if print_all_failures:
                print(count,"(already found up to a constant)")
            continue
        else:
            found_vectors.append(non_constant_part)
        
        # Find actual vector, as an IntegerMatrix object.
        vector = IntegerMatrix(1,degree+compression)
        for coord in range(degree+1):
            vector[0].addmul(A[coord],int(j[coord]))
            
        # Turn into a Python list
        new_vector = [0]*(degree+compression)
        for coord in range(degree+compression):
            new_vector[coord] = vector[0,coord]
        
        if max(new_vector)-min(new_vector)+1 > degree+compression:
            if show_only_dyn_comp:
                print(count,"(doesn't dynamically compress)")
                continue
        else:
            found_dyn_comp = True
        
        coordinates = tuple(map(int,j))
        print(count,"Coordinates:",coordinates,"\tCompression:",degree+compression,"to",max(new_vector)-min(new_vector)+1,end="")
        if show_vector:
            print("\tVector:",new_vector,end="")
        if show_polynomial:
            print("\tPolynomial:",polynomial_tools.poly_coords_to_string(coordinates),end="")
        if show_translated_poly:
            y_trans = 1-min(new_vector)
            print("\tTranslated polynomial:",polynomial_tools.poly_coords_to_string(coordinates,1,y_trans),end="")
        print()
    print("----------")
    return found_dyn_comp

print_dynamical_compression(2,6,100000,show_vector=False,show_polynomial=True,show_translated_poly=True)