from fpylll import *
from fpylll.fplll.gso import MatGSO
import numpy as np
import math

def LambdaDEBasis(d,e):
    basis = IntegerMatrix(d+1,d+e+1)
    for i in range(d+1):
        for j in range(d+e+1):
            basis[i,j] = math.comb(j, i)
    return basis
def PossibleDynCompressing(degree,compression,number_of_solutions):
    # Input: degree, compression, #of solutions
    # Output: coordinates of each polynomial which can possibly send [degree+compression] to [degree+compression]. note not all of these vectors will work!
    if compression < 1:
        print("Invalid compression!")
        return
    
    A = LambdaDEBasis(degree,compression-1)
    M = MatGSO(A)
    _ = M.update_gso()
    enum = Enumeration(M,strategy = EvaluatorStrategy.BEST_N_SOLUTIONS,nr_solutions = number_of_solutions,sub_solutions=False)
    #print((degree+compression)*floor((degree+compression)/2)^2)
    e1 = enum.enumerate(0, degree+1, (degree+compression)*((math.floor((degree+compression)/2))**2), 0)
    # print((math.floor((degree+compression)/2))**2)
    count = 0
    found_vectors = []
    for i,j in e1:
        count += 1
        # Check for unwanted cases
        if max(j[2:]) == 0 and min(j[2:]) == 0:
            # print(count,"(linear solution)")
            continue
        if j[-1] == 0:
            print(count,"(lower degree solution)")
            continue
        #print((i),tuple(map(int,j)))
        
        # Check if we found the same polynomial up to a constant
        non_constant_part = j[1:]
        if non_constant_part in found_vectors:
            # print(count,"(already found up to a constant)")
            continue
        else:
            found_vectors.append(non_constant_part)
        
        # Find actual vector
        vector = IntegerMatrix(1,degree+compression)
        for coord in range(degree+1):
            vector[0].addmul(A[coord],int(j[coord]))
            #print(A[coord],int(j[coord]))
            #print(vector)
            
        # Turn into a Python list
        new_vector = [0]*(degree+compression)
        for coord in range(degree+compression):
            new_vector[coord] = vector[0,coord]
        #print(new_vector)
        
        if max(new_vector)-min(new_vector)+1 <= degree+compression:
            print(poly(tuple(map(int,j)),new_vector),"\tCompression:",degree+compression,"to",max(new_vector)-min(new_vector)+1)


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


def print_poly(poly_v):         # returns the polynomial with coefficients given by poly_v into a "nice" format, as a string
    d = len(poly_v)-1
    
    # poly_v = poly_v*math.factorial(d)  # multiply by d! so that the coeff are integers
    gcd = GCD(poly_v)
    poly_v = poly_v/gcd         # simplifies the polynomial
    poly_string = "("
    for i in range(d+1):        # prints term-by-term, in our nice format
        if poly_v[i]:
            if i==0:
                poly_string+=str(int(poly_v[i])) # prints out the constant term
            elif i==1:
                if poly_v[i] == 1:
                    poly_string+=str("+x") # prints out the x term without the power and without the coeff, but with a plus sign if it is 1
                elif poly_v[i] == -1:
                    poly_string+=str("-x") # prints out the x term without the power and without the coeff, but with a minus sign if it is -1
                else:
                    poly_string+=str(('{:+}'.format(int(poly_v[i])))+"x") # prints out the x term without the power
            else:
                if poly_v[i] == 1:
                    poly_string+=str("+x^"+str(i)) # prints out the x^i term without the coeff, but with a plus sign if it is 1
                elif poly_v[i] == -1:
                    poly_string+=str("-x^"+str(i)) # prints out the x^i term without the coeff, but with a minus sign if it is -1
                else:
                    poly_string+=str(('{:+}'.format(int(poly_v[i])))+"x^"+str(i)) # prints out the x^i term
    poly_string+=str(")/"+str(int(math.factorial(d)/gcd)))          # finish it by dividing by d!, so that the expression looks nicer
    return poly_string

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

    # print("Untranslated polynomial is: ",end="")
    # print_poly(poly_v)

    min_val = min(vals)
    compress_poly = translate(poly_v,min_val)       # translates the polynomial above to get the actual dynamically compressing polynomial
    # print("This translates to dynamically compressing polynomial: ",end="")
    return print_poly(compress_poly)
    # print("          --------              ")


print("Hello World!") # So you know the program is running
for e in range(1,7):
    PossibleDynCompressing(2,e,2000)
