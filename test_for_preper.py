from fractions import Fraction
import sympy
import math

c = 0

# Find a theoretical escape radius for the polynomial z^2 + c
# This is quarantined here in case we have a better way, so we don't have to look for every instance of this.
def find_escape_radius(c):
    return (1 + math.sqrt(1+4*abs(c)))/2

# Find the possible denominators for the polynomial z^2 + c
# This is quarantined here in case we have a better way, so we don't have to look for every instance of this.
def possible_preper_denominators(c):
    if c == 0: return [1]
    if is_square(c.denominator):
        return sympy.divisors(math.isqrt(c.denominator))
    else:
        #print("Denominator is not a square!")
        return
    
# Decides if a non-negative integer x is a perfect square.
def is_square(x):
    if x >= 0:
        return x == math.isqrt(x)**2
    else:
        return False


# Returns potentially preperiodic points, based on escape radius and possible denominators.
def find_potential_preper(c):
    escape_radius = find_escape_radius(c)
    possible_denoms = possible_preper_denominators(c)
    possible_preper = []
    if possible_denoms == None:
        return
    for denominator in possible_denoms:
        numerator_bound = int(escape_radius*denominator) + 1 #to mitigate rounding error
        for numerator in range(-numerator_bound,numerator_bound+1):
            if math.gcd(numerator,denominator) == 1:
                possible_preper.append(Fraction(numerator,denominator))
    possible_preper.sort()
    return possible_preper

# Tests the number z for preperiodicity in z**2 + c
def test_for_preperiodicity(c,z):
    achieved_values = []
    #print(z)
    escape_radius = find_escape_radius(c)
    possible_denoms = possible_preper_denominators(c)
    while z not in achieved_values:
        achieved_values.append(z)
        z = z**2 + c
        #print(z)
        if abs(z) > escape_radius:
            #print("No")
            return False
        if z.denominator not in possible_denoms:
            #print("No")
            return False
    #print("Yes")
    return True

def find_preperiodic_points(c):
    potential_preper = find_potential_preper(c)
    confirmed_preper = []
    for z in potential_preper:
        if test_for_preperiodicity(c,z):
            confirmed_preper.append(z)
    return confirmed_preper

def scan_possible_c(range_min,range_max,denom_min,denom_max):
    for denom in range(denom_min,denom_max+1):
        #print("denom:",denom)
        for numer in range(range_min*denom*denom,range_max*denom*denom+1):
            if math.gcd(numer,denom) == 1:
                found_preper = find_preperiodic_points(Fraction(numer,denom*denom))
                if len(found_preper) > 4:
                    print(Fraction(numer,denom*denom),len(found_preper))
    
scan_possible_c(-10,0,1,50)