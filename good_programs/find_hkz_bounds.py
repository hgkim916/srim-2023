from fpylll import *
import math
import copy
import numpy as np

def lambda_d_e_basis(d,e):
    '''
    Returns the basis (u_0, u_1, ..., u_(d+1)) as an IntegerMatrix object, with rows as basis vectors.
    '''
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def hkz_reduction(matrix):
    '''
    Returns an HKZ-reduced basis for the given lattice basis, where rows are basis vectors.
    '''
    # Note: this is not documented properly, but BKZ.reduction changes the initial matrix too.
    matrix_toreduce = copy.copy(matrix) # By making a copy first, the original matrix is not affected.
    return BKZ.reduction(matrix_toreduce, BKZ.Param(block_size = matrix_toreduce.nrows))

def hkz_succ_min_bounds_euclidean(matrix,minima):
    '''
    Computes lower and upper bounds for the Euclidean norm of the given successive minima,
      given a matrix, by HKZ-reducing the matrix.
    The results come from Theorem 2 in Schnorr (1994).
    Keep in mind this function may be subject to floating-point precision errors (which may be avoidable by sticking to integers - sorry!).

    Parameters:
        matrix (IntegerMatrix): Basis for the lattice, where rows are basis vectors.
        minima (tuple): The successive minima for which to compute bounds.

    Returns:
        Two dictionaries: lower_bounds, upper_bounds.
    '''
    matrix_reduced = hkz_reduction(matrix)
    lower_bounds = {}
    upper_bounds = {}
    for i in minima:
        b_i = matrix_reduced[i-1] # Find the appropriate vector
        approximation_factor = math.pow(4/(i+3),0.5)
        lower_bounds[i] = abs(b_i)*approximation_factor
        upper_bounds[i] = abs(b_i)/approximation_factor
    return lower_bounds,upper_bounds

# Modifies the results of HKZ_SuccessiveMinimaBounds_Euclidean to obtain lower and upper bounds for
# the Max norm.
# Note if we have an integer lattice then we can take ceilings and floors.
def hkz_succ_min_bounds_maxnorm(matrix,minima):
    '''
    Computes lower and upper bounds for the max-norm of the given successive minima,
      given a matrix, by HKZ-reducing the matrix, and using Euclidean bounds.
    The results come from Theorem 2 in Schnorr (1994).
    Keep in mind this function may be subject to minor floating-point precision errors (which may be avoidable by sticking to integers - sorry!).

    Parameters:
        matrix (IntegerMatrix): Basis for the lattice, where rows are basis vectors.
        minima (tuple): The successive minima for which to compute bounds.

    Returns:
        Two dictionaries: lower_bounds, upper_bounds.
    '''
    dimension = len(matrix[0])
    lower_bounds,upper_bounds = hkz_succ_min_bounds_euclidean(matrix,minima)
    for i in minima:
        lower_bounds[i] = lower_bounds[i]/math.pow(dimension,0.5)
    return lower_bounds,upper_bounds

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
