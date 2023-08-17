from fpylll import *
import math
import copy
import numpy as np

# The basis for LambdaDE in terms of row vectors.
def LambdaDEBasis(d,e):
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

# HKZ reduction on the given matrix.
def HKZReduction(matrix):
    # Note: this is not documented properly, but BKZ.reduction changes the initial matrix too.
    matrix_toreduce = copy.copy(matrix) # By making a copy first, the original matrix is not affected.
    return BKZ.reduction(matrix_toreduce, BKZ.Param(block_size = matrix_toreduce.nrows))

# Uses Theorem 2 in Schnorr (1994) to return lower and upper bounds for the given successive minima.
# The input variable "minima" can be a list/tuple of positive integers.
# Returns two dicts: lower_bounds, upper_bounds.
def HKZ_SuccessiveMinimaBounds_Euclidean(matrix,minima):
    matrix_reduced = HKZReduction(matrix)
    lower_bounds = {}
    upper_bounds = {}
    for i in minima:
        # Find the appropriate vector
        b_i = matrix_reduced[i-1]
        approximation_factor = math.pow(4/(i+3),0.5)
        lower_bounds[i] = abs(b_i)*approximation_factor
        upper_bounds[i] = abs(b_i)/approximation_factor
    return lower_bounds,upper_bounds

# Modifies the results of HKZ_SuccessiveMinimaBounds_Euclidean to obtain lower and upper bounds for
# the Max norm.
# Note if we have an integer lattice then we can take ceilings and floors.
def HKZ_SuccessiveMinimaBounds_MaxNorm(matrix,minima):
    dimension = len(matrix[0])
    lower_bounds,upper_bounds = HKZ_SuccessiveMinimaBounds_Euclidean(matrix,minima)
    for i in minima:
        lower_bounds[i] = lower_bounds[i]/math.pow(dimension,0.5)
    return lower_bounds,upper_bounds


# Testing/Examples. Can delete.
A = LambdaDEBasis(10,20)
print(A)
A_reduced = HKZReduction(A)
print(A)
print(A_reduced)
print(len(A_reduced[1]))
print(abs(A_reduced[1]))
print(A[1].norm())
B,C = HKZ_SuccessiveMinimaBounds_Euclidean(A,[1,2])
print(B)
print(C)
print(HKZ_SuccessiveMinimaBounds_MaxNorm(A,[1,2,3]))
