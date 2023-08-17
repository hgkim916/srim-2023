from fpylll import *
from fpylll.fplll.gso import MatGSO
import math

def LambdaDEBasis(d,e):
    '''Returns the basis (u_0, u_1, ..., u_(d+1)) as an IntegerMatrix object, with rows as basis vectors.
    '''
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def PossibleDynCompressing(degree,compression,number_of_solutions,
                           show_only_dyn_comp = True,
                           show_vector = True,
                           verbose = True):
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
        verbose (bool): (Default: True) Display more output information. If false, can cut down on the amount of outsput significantly.
    '''
    # First check for validity of input
    if compression < 1:
        raise ValueError("Invalid compression!")
    
    print("----------PossibleDynCompressing")
    print("Looking for IVPs of degree",degree,"which compress ["+str(degree+compression)+"]")

    # This part is basically taken from fpylll's Github page. We're not exactly sure how it works.
    A = LambdaDEBasis(degree,compression-1)
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
            if verbose:
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
        
        if show_vector:
            print(count,"Coordinates:",tuple(map(int,j)),"\tCompression:",degree+compression,"to",max(new_vector)-min(new_vector)+1,"\tVector:",new_vector)
        else:
            print(count,"Coordinates:",tuple(map(int,j)),"\tCompression:",degree+compression,"to",max(new_vector)-min(new_vector)+1)
    print("----------")
    return found_dyn_comp

#PossibleDynCompressing(10,10,10000,show_only_dynamically_compressing = True,show_vector = True)