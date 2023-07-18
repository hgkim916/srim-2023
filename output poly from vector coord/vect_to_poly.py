import numpy as np
import math

def unsigned_stirling_num_recc(n,k):        # computes the unsigned Stirling numbers of the first kind recursively, see wiki page below
    if k == 0 and n == 0:
        return 1
    elif k == 0 or n == 0:
        return 0
    else:
        return (n-1)*unsigned_stirling_num_recc(n-1,k) + unsigned_stirling_num_recc (n-1,k-1)


def basis_poly(i,d):      # returns the polynomial (x choose i) as a vector of the form [a_0, a_1, ..., a_{d-1}, a_d]
    # Uses the expansion provided by the following wikipedia page: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
    p = np.zeros(d+1)
    for j in range(i+1):
        p[j] = (-1)**(i-j)*unsigned_stirling_num_recc(i,j)/math.factorial(i)
    return p

def output_poly(v):             # outputs the polynomial given by the vector v
    d = len(v)-1
    poly_v = np.zeros(d+1)
    for i in range(len(v)):
        poly_v += v[i]*basis_poly(i,d)      # stores the polynomial coefficients in the form [a_0, a_1, ..., a_{d-1}, a_d]

    poly_v = poly_v*math.factorial(d)  # multiply by d! so that the coeff are integers
    print("(",end="")
    for i in range(d+1):
        if poly_v[i]:
            print(('{:+}'.format(int(poly_v[i])))+"x^"+str(i),end='') # prints out the x^i term
    print(")/"+str(math.factorial(d)))          # finish it by dividing by d!, so that the expression looks nicer



# user inputs
v = [3,-3,1]            # the vector we are trying to simplify, with coordinates w.r.t the basis u_0, u_1 etc of the lattice
output_poly(v)