#!/usr/bin/env python3
# ^ the above code makes it run properly on the server (I think?)

from fpylll import LLL,BKZ, IntegerMatrix  # Lattice reduction algorithms
import math         # For square rooting and stuff like that
import numpy as np  # For various matrix-related things
from matplotlib import pyplot as plt # For the plots

# WHAT THIS PROGRAM DOES:
# For all the values of d in the given range, it looks at e in [0,40] and saves into a file 
# the last value of e for which the third basis vector of the HKZ reduced basis is not 
# equivalent to the optimal quadratic (upon translation and reflection)


def LambdaDEBasis(d,e):     # This creates the lattice \Lambda_{d,e}
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def LLL_reduced_b3(d,e):           # This returns the 3rd basis vector of the LLL reduced basis of the lattice \Lambda_{d,e}
    A = LambdaDEBasis(d,e)
    W = LLL.Wrapper(A)      
    W()                         # Runs the reduction
    if not LLL.is_reduced(A):   # This is just a sanity check, can be removed
        print("WARNING! Sanity check failed, lattice is not LLL reduced")
    # print(A)
    return A[2]                    # We don't need the other basis vectors

def HKZ_reduced_b3(d,e):           # This returns the 3rd basis vector of the HKZ reduced basis of the lattice \Lambda_{d,e}
    A = LambdaDEBasis(d,e)
    A_red = BKZ.reduction(A, BKZ.Param(d+1))    # Runs BKZ reduction with beta = d+1, which is the same as HKZ reduction
    # print(A_red)
    return A_red[2]                # We don't need the other basis vectors

def q_min(x):                      # This is the optimal quadratic (up to translation and reflection)
    return (-x**2+x)/2



# print("Hello World!") # So you know the program is running

# Set details of the range of d we want to investigate:
minsize = 2
maxsize = 70
step = 1

# We will output the results to a file in the subfolder Results
f = open("./Results/is_b3_quad_for_large_e.txt","w") # We open and close the file to delete any previous results
f.close()

for d in range(minsize,maxsize+1):

    max_e_not_agrees = -1       # We will store in this the maximum value of e for which the reduced b_3 does NOT give the quadratic q_min
    es = range(41)          # we investigate 0<=e<=40
    for i in range(len(es)):
        b3_HKZ = HKZ_reduced_b3(d,es[i])
        y_translated = b3_HKZ - min(b3_HKZ)*np.ones(len(b3_HKZ))    # translate the values so that they are non-negative
        translation_in_x = np.argmin(y_translated)          # by how much we need to shif in x to get the minimum between 0 and 1 i.e. to match (-x^2+x)/2
        for j in range(len(b3_HKZ)):
            if y_translated[j] != q_min(j - translation_in_x) and y_translated[j] != -q_min(j - translation_in_x):
                max_e_not_agrees = es[i]      # the quadratic doesn't agree with our b_3, and hence we increase the maximum
                break   # no point checking further, exit and move on to next value of e
                   
    f = open("./Results/is_b3_quad_for_large_e.txt","a")     # We now open the file in append mode, to add the results to the end        
    f.write("# d = "+str(d)+'\n')                            # Write results to file
    f.write("Largest value of e that doesn't have b_3 as a translation of the optimal quadratic is e = "+str(max_e_not_agrees)+"\n") 
    f.close()                                                # Close the file