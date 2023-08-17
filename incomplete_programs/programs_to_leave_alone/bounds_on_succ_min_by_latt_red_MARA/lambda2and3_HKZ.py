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

def LLL_reduced(d,e):           # This returns the LLL reduced basis (as row vecs) of the lattice \Lambda_{d,e}
    A = LambdaDEBasis(d,e)
    W = LLL.Wrapper(A)      
    W()                         # Runs the reduction
    if not LLL.is_reduced(A):   # This is just a sanity check, can be removed
        print("WARNING! Sanity check failed, lattice is not LLL reduced")
    # print(A)
    return A                    # Returns the entire basis. Could be adjusted to only return e.g. a specific row vector

def HKZ_reduced(d,e):           # This returns the HKZ reduced basis (as row vecs) of the lattice \Lambda_{d,e}
    A = LambdaDEBasis(d,e)
    A_red = BKZ.reduction(A, BKZ.Param(d+1))    # Runs BKZ reduction with beta = d+1, which is the same as HKZ reduction
    # print(A_red)
    return A_red                # Returns the entire basis. Could be adjusted to only return e.g. a specific row vector

print("Hello World!") # So you know the program is running

# Set details of the range of d we want to investigate:
minsize = 10           
maxsize = 65
step = 5

for d in range (minsize,maxsize+1,step):

    print(d)        # See progress the program is making

    es = range(math.floor(math.log2(d/2))+10)       # plotting for e from 0 to what they proved in the paper and another 9 values after that

    maxReqBound3 = [(d+e+1)/2 for e in es]          # we want lambda3 to be at most this for dynamical compression to even exist
    # maxReqBound3 = [(d+e)/2 for e in es]            # what I expect the bound we want on lambda2 to be, altough I haven't verified it
    
    lhsB2 = np.zeros(len(es))   
    rhsB2 = np.zeros(len(es))                  
    lhsB3 = np.zeros(len(es))
    rhsB3 = np.zeros(len(es))               # Initialise the vectors that we will plot, so we can access without worry afterwards

    for i in range(len(es)):
        latt = HKZ_reduced(d,es[i])         # get reduced basis

        # Calculate lower and upper bounds, as given in the Overleaf document:
        lhsB2[i] = math.sqrt((4)/(5*(d+es[i]+1))) * np.linalg.norm(latt[1],2)       # latt[1] is b_2
        rhsB2[i] =  math.sqrt((5)/(4)) * np.linalg.norm(latt[1],2)
        lhsB3[i] = math.sqrt((2)/(3*(d+es[i]+1))) * np.linalg.norm(latt[2],2)       # latt[2] is b_3
        rhsB3[i] =  math.sqrt((3)/(2)) * np.linalg.norm(latt[2],2)  

    # Plot the graph - should be fairly self-explanatory
    plt.plot(es,maxReqBound3,'y--',label="max value lambda_3 can take")
    plt.plot(es,lhsB2,'g',label="lower bound given by b_2")
    plt.plot(es,rhsB2,'g',label="upper bound given by b_2")
    plt.plot(es,lhsB3,'k',label="lower bound given by b_3")
    plt.plot(es,rhsB3,'k',label="upper bound given by b_3")
    plt.axvline(x = math.floor(math.log2(d/2)), color = 'y',label="value of e the result has been proven to")
    plt.fill_between(es,lhsB2,rhsB2, alpha=0.2)         # Fill in the possible range for \lambda_2
    plt.fill_between(es,lhsB3,rhsB3, alpha=0.2)         # Fill in the possible range for \lambda_3
    plt.legend(loc="upper left")
    plt.ylim([0, 3*max(maxReqBound3)])
    plt.title("d="+str(d))
    plt.xlabel("e")
    plt.ylabel("value of second and 3rd successive minima")
    plt.savefig("Bounding second and third successive minima by HKZ for d="+str(d))     # Automatically save plot to file
    # plt.show()              # To display the picture, if you're running it in a Jupiter notebook or smth - not useful otherwise
    plt.close()         # Clears the plot