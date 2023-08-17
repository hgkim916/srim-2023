#!/usr/bin/env python3
# ^ the above code makes it run properly on the server (I think?)

from fpylll import LLL,BKZ, IntegerMatrix  # Lattice reduction algorithms
import math         # For square rooting and stuff like that
import numpy as np  # For various matrix-related things
from matplotlib import pyplot as plt # For the plots

# WHAT THIS PROGRAM DOES:
# Given some values of d to look at, it will plot a graph for each one, with values of e on teh x-axis, 
# and values of the bounds we get for the second and third successive minima respectively on the y-axis.
# More about the bounds can be found on the Overleaf write-up titled 
#                       "Using lattice reduction algorithms to bound the successive minima"


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

print("Hello World!") # So you know the program is running

# Set details of the range of d we want to investigate:
minsize = 60
maxsize = 70
step = 1

for d in range (minsize,maxsize+1,step):

    print(d)        # See progress the program is making

    es = range(math.floor(math.log2(d/2))+10)       # plotting for e from 0 to what they proved in the paper and another 9 values after that

    maxReqBound = [(d+e+1)/2 for e in es]          # we want lambda3 to be at most this for dynamical compression to even exist

    delta = 0.99                                #
    alpha = 1/(delta-0.25)                      # Constants of the LLL reduction, used in our bound computation

    lhsLLLred = np.zeros(len(es))               #
    rhsLLLred = np.zeros(len(es))               # 
    lhsHKZred = np.zeros(len(es))               # Initialise the vectors that we will plot, so we can access without worry afterwards
    rhsHKZred = np.zeros(len(es))               #

    for i in range(len(es)):

        b3_LLL = LLL_reduced_b3(d,es[i])        # Get the third LLL reduced basis vector

        # Calculate lower and upper LLL bounds, as given in the Overleaf document:
        lhsLLLred[i] = math.sqrt((math.pow(alpha,-d))/(d+es[i]+1)) * np.linalg.norm(b3_LLL,2)
        rhsLLLred[i] = alpha * np.linalg.norm(b3_LLL,2)

        b3_HKZ = HKZ_reduced_b3(d,es[i])        # Get the third HKZ reduced basis vector

        # Calculate lower and upper HKZ bounds, as given in the Overleaf document:
        lhsHKZred[i] = math.sqrt((2)/(3*(d+es[i]+1))) * np.linalg.norm(b3_HKZ,2)
        rhsHKZred[i] =  math.sqrt((3)/(2)) * np.linalg.norm(b3_HKZ,2)

    # Plot the graph - should be fairly self-explanatory
    plt.plot(es,maxReqBound,'g--',label="max value lambda_3 can take")
    plt.plot(es,lhsLLLred,'k',label="lower bound given by LLL")
    plt.plot(es,rhsLLLred,'k',label="upper bound given by LLL")
    plt.plot(es,lhsHKZred,'b',label="lower bound given by HKZ")
    plt.plot(es,rhsHKZred,'b',label="upper bound given by HKZ")
    plt.axvline(x = math.floor(math.log2(d/2)), color = 'y',label="value of e the result has been proven to")
    plt.fill_between(es,lhsHKZred,rhsHKZred, alpha=0.2)# Fill in the possible range for \lambda_3 - only do for HKZ as it's just better
    plt.legend(loc="upper left")
    plt.ylim([0, 3*max(maxReqBound)])
    plt.savefig("Bounding 3rd succ min for d="+str(d))              # Automatically save plot to file
    # plt.show()              # To display the picture, if you're running it in a Jupiter notebook or smth - not useful otherwise
    plt.close()             # Clears the plot