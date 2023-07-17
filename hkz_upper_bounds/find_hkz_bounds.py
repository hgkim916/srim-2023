#!/usr/bin/env python3

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


def HKZ_Bounds_LambdaDE(d,e,minimum):
    A = LambdaDEBasis(d,e)
    bounds = HKZ_SuccessiveMinimaBounds_MaxNorm(A,[minimum])
    return bounds[0][minimum],bounds[1][minimum]

def Print_Bounds_d_e(d,max_e):
    for e in range(1,max_e):
        print(d,e,HKZ_Bounds_LambdaDE(d,e,3)[0],math.floor((d+e+1)/2))

def Find_emax_UpperBound(d, print_progress = False, print_result = True):
    max_e = 10
    e = 1
    while e <= max_e:
        lower_bound = HKZ_Bounds_LambdaDE(d,e,3)[0]
        if print_progress:
            print(d,e,lower_bound,math.floor((d+e+1)/2),lower_bound<math.floor((d+e+1)/2))
        if lower_bound<math.floor((d+e+1)/2):
            max_e = e + 4 # Require 4 Falses in a row to give up
        e += 1
    if print_result:
        print("----------------",d,max_e - 4)
    return max_e - 4

Find_emax_UpperBound(20,print_progress = True)
