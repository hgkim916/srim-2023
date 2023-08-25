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