# This program implements an enumeration algorithm for finding all vectors in a given Euclidean ball.
# However, it doesn't work very reliably because of round-off errors.

import numpy
import math

# Initilizes a basis for LambdaDE into an array.
# The basis vectors in this case are rows of the matrix.
def LambdaDEBasis_Array(d,e):
    matrix = numpy.zeros((d+1,d+e+1))
    for i in range(d+1):
        for j in range(d+e+1):
            matrix[i][j] = math.comb(j,i)
    return matrix

##B = LambdaDEBasis_Array(2,3)
##A = numpy.matmul(B,numpy.transpose(B))
##R = numpy.transpose(numpy.linalg.cholesky(A))
##Rinv = numpy.linalg.inv(R)
##print(B)
##print(A)
##print(R)
##print(Rinv)

# Calculates every vector shorter than C in the lattice with basis of row vectors B = basis.
# This is Algorithm 2.8 in Fincke/Pohst 1985. This is space efficient but very time-inefficient.
# The biggest hurdle with this algorithm is that there is a lot of numerical error.
def Enumerate_Naive(B,C):
    m = len(B)
    A = numpy.matmul(B,numpy.transpose(B))
    R = numpy.transpose(numpy.linalg.cholesky(A))
    Q = numpy.zeros((m,m))
    for i in range(m):
        Q[i][i] = (R[i][i])**2
        for j in range(i+1,m):
            Q[i][j] = R[i][j]/R[i][i]
    print(Q)
    T = numpy.zeros(m)
    U = numpy.zeros(m)
    x = numpy.zeros(m)
    UBx = numpy.zeros(m)
    i = m-1
    T[i] = C
    U[i] = 0
    keep_going = True
    count = 0
    while keep_going: # Step 2
##        print("Starting step 2")
##        print(x)
##        print(T)
##        print(U)
##        print(i)
        Z = (T[i]/Q[i][i])**0.5
        UBx[i] = math.floor(Z - U[i])
        x[i] = math.ceil(-Z - U[i]) -1
        while True: # Step 3
##            print("Starting step 3")
            x[i]+=1
            if x[i] > UBx[i]:
##                print("Starting step 4")
                i+=1 # Step 4
            elif i != 0: # Step 5
##                print("Starting step 5")
                i-=1
                U[i] = 0
                for j in range(i+1,m):
                    U[i] += Q[i][j]*x[j]
                T[i] = T[i+1] - Q[i+1][i+1]*(x[i+1] + U[i+1])**2
                break
            else: # Step 6
##                print("Starting step 6")
                if max(x) == 0 and min(x) == 0:
                    keep_going = False
                    break
                else:
                    count += 1
                    print(count,x,"Success!?",int(C - T[0] + Q[0][0]*(x[0]+U[0])**2))
    


B = LambdaDEBasis_Array(3,7)
print(B)
A = numpy.matmul(B,numpy.transpose(B))
print(A)
A = Enumerate_Naive(B,1000)

