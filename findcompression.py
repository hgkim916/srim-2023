# To run: go to https://sagecell.sagemath.org/ or download SageMath (lmao why is it so hard to)
from fpylll import *
from fpylll.fplll.gso import MatGSO
import math

def LambdaDEBasis(d,e):
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def PossibleDynCompressing(degree,compression,number_of_solutions,show_only_dynamically_compressing = False,show_vector = False):
    # Input: degree, compression, #of solutions
    # Output: coordinates of each polynomial which can possibly send [degree+compression] to [degree+compression]. note not all of these vectors will work!
    if compression < 1:
        print("Invalid compression!")
        return
    
    A = LambdaDEBasis(degree,compression-1)
    M = MatGSO(A)
    _ = M.update_gso()
    enum = Enumeration(M,strategy = EvaluatorStrategy.BEST_N_SOLUTIONS,nr_solutions = number_of_solutions,sub_solutions=False)
    #print((degree+compression)*floor((degree+compression)/2)^2)
    e1 = enum.enumerate(0, degree+1, (degree+compression)*(math.floor((degree+compression)/2)**2), 0)
    count = 0
    found_vectors = []
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
            print(count,"(already found up to a constant)")
            continue
        else:
            found_vectors.append(non_constant_part)
        
        # Find actual vector
        vector = IntegerMatrix(1,degree+compression)
        for coord in range(degree+1):
            vector[0].addmul(A[coord],int(j[coord]))
            #print(A[coord],int(j[coord]))
            #print(vector)
            
        # Turn into a Python list
        new_vector = [0]*(degree+compression)
        for coord in range(degree+compression):
            new_vector[coord] = vector[0,coord]
        #print(new_vector)
        
        if show_only_dynamically_compressing and max(new_vector)-min(new_vector)+1 > degree+compression:
            print(count,"(doesn't dynamically compress)")
            continue
        
        print(count,"Coordinates:",tuple(map(int,j)),"\tCompression:",degree+compression,"to",max(new_vector)-min(new_vector)+1,"\tVector:",new_vector)
        
print("Hello World!") # So you know the program is running
PossibleDynCompressing(10,10,10000,show_only_dynamically_compressing = True,show_vector = True)

