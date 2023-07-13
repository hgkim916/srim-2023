#!/usr/bin/env python3
from fpylll import LLL,BKZ, IntegerMatrix, GSO
import math
import fpylll
import numpy as np
from matplotlib import pyplot as plt

def LambdaDEBasis(d,e):
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis

def LLL_reduced(d,e):
    A = LambdaDEBasis(d,e)
    W = LLL.Wrapper(A)
    W()
    if not LLL.is_reduced(A):
        print("WARNING! Sanity check failed, lattice is not LLL reduced")
    # print(A)
    return A

def HKZ_reduced(d,e):
    A = LambdaDEBasis(d,e)
    A_red = BKZ.reduction(A, BKZ.Param(d+1))
    # if not LLL.is_reduced(A):
    #     print("WARNING! Sanity check failed, lattice is not LLL reduced")
    print(A_red)
    return A_red

print("Hello World!") # So you know the program is running
minsize = 10
maxsize = 65
step = 5
for d in range (minsize,maxsize+1,step):
    print(d)
    es = range(math.floor(math.log2(d/2))+10)
    maxReqBound2 = [(d+e+1)/2 for e in es]
    # maxReqBound3 = [(d+e)/2 for e in es]
    lhsB2 = np.zeros(len(es))
    rhsB2 = np.zeros(len(es))
    lhsB3 = np.zeros(len(es))
    rhsB3 = np.zeros(len(es))
    for i in range(len(es)):
        latt = HKZ_reduced(d,es[i])
        lhsB2[i] = math.sqrt((4)/(5*(d+es[i]+1))) * np.linalg.norm(latt[1],2)
        rhsB2[i] =  math.sqrt((5)/(4)) * np.linalg.norm(latt[1],2)
        lhsB3[i] = math.sqrt((2)/(3*(d+es[i]+1))) * np.linalg.norm(latt[2],2)
        rhsB3[i] =  math.sqrt((3)/(2)) * np.linalg.norm(latt[2],2)  


    plt.plot(es,maxReqBound2,'y--',label="max value lambda_3 can take")
    plt.plot(es,lhsB2,'g',label="lower bound given by b_2")
    plt.plot(es,rhsB2,'g',label="upper bound given by b_2")
    plt.plot(es,lhsB3,'k',label="lower bound given by b_3")
    plt.plot(es,rhsB3,'k',label="upper bound given by b_3")
    plt.axvline(x = math.floor(math.log2(d/2)), color = 'y',label="value of e the result has been proven to")
    plt.fill_between(es,lhsB2,rhsB2, alpha=0.2)
    plt.fill_between(es,lhsB3,rhsB3, alpha=0.2)
    plt.legend(loc="upper left")
    plt.ylim([0, 3*max(maxReqBound2)])
    plt.title("d="+str(d))
    plt.xlabel("e")
    plt.ylabel("value of second and 3rd successive minima")
    plt.savefig("Bounding second and third successive minima by HKZ for d="+str(d))
    plt.close()