import numpy as np
import math

def unsigned_stirling_first_kind(n,k):
    '''
    Computes the unsigned Stirling numbers of the first kind, recursively.
    See the Wikipedia page for Stirling numbers of the 1st kind for more information.
    '''
    if k == 0 and n == 0:
        return 1
    elif k == 0 or n == 0:
        return 0
    else:
        return (n-1)*unsigned_stirling_first_kind(n-1,k) + unsigned_stirling_first_kind (n-1,k-1)

def basis_poly_coeffs(i,d):
    '''
    Returns the coefficients of the polynomial d!*(x choose i), as a list [a_0, a_1, ..., a_{d-1}, a_d], where a_j is the coefficient of x^j.
    See the Wikipedia page for Stirling numbers of the 1st kind for more information.

    Requires that d is greater than or equal to i.
    '''
    if d < i:
        raise ValueError("d must be greater than or equal to i.")

    p = np.zeros(d+1,dtype=np.int_)
    for j in range(i+1):
        p[j] = (-1)**(i-j)*unsigned_stirling_first_kind(i,j)*math.perm(d,d-i)
    return p

def poly_coeffs_to_string(coefficients,denom = 1):
    '''
    Given a list of coefficients, poly_v = d!*[a_0, a_1, ..., a_d], returns the polynomial
    (a_0+a_1x+...+a_dx^d)/denom, simplified, as a string.
    '''
    d = len(coefficients)-1

    poly_string = "("
    for i in range(d+1):        # go term-by-term
        if coefficients[i] == 0: continue # skip term if 0

        if i==0:
            poly_string += str(coefficients[i])
            continue

        if coefficients[i] == 1:
            poly_string+="+x"
        elif coefficients[i] == -1:
            poly_string+="-x"
        else:
            poly_string+='{:+}'.format(coefficients[i])+"x"

        if i == 1:
            continue
        else:
            poly_string+="^"+str(i)
    
    if denom == 1:
        poly_string+=")"
    else:
        poly_string+=")/"+str(denom)
    
    return poly_string

def translated_poly_coeffs(coefficients,x_translation,y_translation):
    '''
    Given the coefficients of a polynomial p, returns the coefficients of p(x-x_translation)+y_translation.
    Note the sign of x_translation! We're translating positively to the right and upwards.
    '''
    if x_translation == 0:
        coefficients[0] += y_translation
        return coefficients
    
    d = len(coefficients)-1
    translated_coeffs = np.zeros(d+1,dtype=np.int_)
    
    for j in range(d+1):
        for i in range(j,d+1):
            if x_translation > 0: # give p(x-1)
                translated_coeffs[j]+=((-1)**(i-j))*coefficients[i]*math.comb(i,j)       # the x^j term in a_i*(x-1)^i  (where p(x) = a_0 + ... + a_d*x^d)
            else:
                translated_coeffs[j]+=coefficients[i]*math.comb(i,j)

    if abs(x_translation) == 1:
        translated_coeffs[0] += y_translation
        return translated_coeffs 
    elif x_translation > 1:
        return translated_poly_coeffs(translated_coeffs,x_translation-1,y_translation)
    else:
        return translated_poly_coeffs(translated_coeffs,x_translation+1,y_translation)

def poly_coords_to_coeffs(coords):
    '''
    Given coords [u_0,u_1,...,u_d], gives the coefficients [a_0,a_1,...,a_d]
    of the polynomial u_0(x choose 0)+u_1(x choose 1)+...+u_d(x choose d), where a_i is the coefficient of x^i.
    '''
    d = len(coords)-1
    coefficients = np.zeros(d+1,dtype=np.int_)
    for i in range(len(coords)):
        coefficients += coords[i]*basis_poly_coeffs(i,d)      # stores the polynomial coefficients in the form d!*[a_0, a_1, ..., a_{d-1}, a_d] 
    return coefficients

def poly_coords_to_string(coords,x_translation=0,y_translation=0):
    '''
    Given coords [u_0,u_1,...,u_d], prints the polynomial u_0(x choose 0)+u_1(x choose 1)+...+u_d(x choose d),
    simplifying to form a_0+a_1x+...+a_dx^d.
    Optionally, also translate the polynomial first.
    '''
    d = len(coords)-1
    coefficients = poly_coords_to_coeffs(coords)
    # Note we need to translate the polynomial in the y-direction by a multiple of d!.
    translated_coeffs = translated_poly_coeffs(coefficients,x_translation=x_translation,y_translation=y_translation*math.factorial(d))
    return poly_coeffs_to_string(translated_coeffs,denom=math.factorial(d))