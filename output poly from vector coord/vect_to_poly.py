# SET-UP - only need to run once

import numpy as np
import math

def unsigned_stirling_num_recc(n,k):        # computes the unsigned Stirling numbers of the first kind recursively, see wiki page below
    if k == 0 and n == 0:
        return 1
    elif k == 0 or n == 0:
        return 0
    else:
        return (n-1)*unsigned_stirling_num_recc(n-1,k) + unsigned_stirling_num_recc (n-1,k-1)

def basis_poly(i,d):      # returns the polynomial d!*(x choose i) as a vector of the form [a_0, a_1, ..., a_{d-1}, a_d]
    # Uses the expansion provided by the following wikipedia page: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
    p = np.zeros(d+1)
    for j in range(i+1):
        p[j] = ((-1)**(i-j)*unsigned_stirling_num_recc(i,j)/math.factorial(i))*math.factorial(d)
    return p

def GCD(v):     # recursively calculates the GCD of all the elements of v
    if len(v)==1:
        return int(v[0])
    else:
        return int(math.gcd(int(v[0]),GCD(v[1:])))

def print_poly(poly_v):         # prints the polynomial with coefficients given by poly_v into a "nice" format
    d = len(poly_v)-1
    
    # poly_v = poly_v*math.factorial(d)  # multiply by d! so that the coeff are integers
    gcd = GCD(poly_v)
    poly_v = poly_v/gcd         # simplifies the polynomial
    print("(",end="")
    for i in range(d+1):        # prints term-by-term, in our nice format
        if poly_v[i]:
            if i==0:
                print(int(poly_v[i]),end='') # prints out the constant term
            elif i==1:
                if poly_v[i] == 1:
                    print("+x",end='') # prints out the x term without the power and without the coeff, but with a plus sign if it is 1
                elif poly_v[i] == -1:
                    print("-x",end='') # prints out the x term without the power and without the coeff, but with a minus sign if it is -1
                else:
                    print(('{:+}'.format(int(poly_v[i])))+"x",end='') # prints out the x term without the power
            else:
                if poly_v[i] == 1:
                    print("+x^"+str(i),end='') # prints out the x^i term without the coeff, but with a plus sign if it is 1
                elif poly_v[i] == -1:
                    print("-x^"+str(i),end='') # prints out the x^i term without the coeff, but with a minus sign if it is -1
                else:
                    print(('{:+}'.format(int(poly_v[i])))+"x^"+str(i),end='') # prints out the x^i term
    print(")/"+str(int(math.factorial(d)/gcd)))          # finish it by dividing by d!, so that the expression looks nicer

def translate(v,min_val):          
    # gives the coefficients of p(x-1)+min_val+1 in terms of x, being given the coefficients of p and the minimum integer value it takes on {0,...,d+e}
    trans_v = np.zeros(len(v))
    d = len(v)-1
    for j in range(len(trans_v)):
        for i in range(j,len(trans_v)):
            trans_v[j]+=((-1)**(i-j))*v[i]*math.comb(i,j)       # the x^j term in a_i*(x-1)^j  (where p(x) = a_0 + ... + a_d*x^d)
    trans_v[0]+=(1-min_val)*math.factorial(d)        # makes sure the min value achieved is precisely 1 by subtractinf the initial min and adding 1
    return trans_v

def poly(v,vals):               # outputs the polynomial given by the vector v
                                # vals is the values taken by this poly on {0,...,d+e}, used in rescaling it
    d = len(v)-1
    poly_v = np.zeros(d+1)
    for i in range(len(v)):
        poly_v += v[i]*basis_poly(i,d)      # stores the polynomial coefficients in the form d!*[a_0, a_1, ..., a_{d-1}, a_d] 
                                            #       (the d! is to prevent rounding errors)

    print("Untranslated polynomial is: ",end="")
    print_poly(poly_v)

    min_val = min(vals)
    compress_poly = translate(poly_v,min_val)       # translates the polynomial above to get the actual dynamically compressing polynomial
    print("This translates to dynamically compressing polynomial: ",end="")
    print_poly(compress_poly)
    print("          --------              ")


# user inputs
v = [3,-3,1]            # the vector we are trying to simplify, with coordinates w.r.t the basis u_0, u_1 etc of the lattice
poly(v)