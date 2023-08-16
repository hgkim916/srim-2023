from fractions import Fraction
import math

def norm(x,y):
    return math.sqrt(x**2 + y**2)

def Henon_quadr_given_ab(a,b,x_coeff = -1):
    def f(x,y):
        return Fraction(y),Fraction(x_coeff*x+(y**2+Fraction(a,b)))
    return f

def get_euclidean_bound_quadr(a,b):     # returns the max Euclidean norm a preperiodic point can have
    return 100*a**2

def get_p_adic_bound_quadr(a,b):        # returns the "worst" denominator a preperiodic point can have (on p-adic grounds)
    return b
    # return 1        # this says we only want integers

def trace_pt(f, eucl_bound,padic_bound,x,y):
    orbit = []
    while [x,y] not in orbit:
        orbit.append([x,y])
        x,y = f(x,y)
        if norm(x,y)>eucl_bound or Fraction(x*padic_bound).denominator>1 or Fraction(y*padic_bound).denominator>1:     
            # Escape, either on Euclidean or padic grounds
            return None
    return orbit


def get_max_cycle_quadr(f,a,b):
    eucl_bound = get_euclidean_bound_quadr(a,b)
    padic_bound = get_p_adic_bound_quadr(a,b)
    max_cycle = 0
    point_on_max = []
    checked_points = []
    # We will search all (x,y) in Q^2 with the bounds above by multiplying by the p_adic bound and searching in a box in Z^2:
    search_space = eucl_bound*padic_bound
    for x in range(-search_space,search_space+1):
        for y in range(-search_space,search_space+1):
            X,Y = Fraction(x,padic_bound),Fraction(y,padic_bound)   # The point we want to trace
            if [X,Y] not in checked_points:         # Do not check if already part of a previously discovered cycle
                orbit = trace_pt(f,eucl_bound,padic_bound,X,Y)
                if orbit != None:
                    checked_points.extend(orbit)        # Add to the list of points we already know about
                    if len(orbit)>max_cycle:
                        max_cycle = len(orbit)
                        point_on_max = [X,Y]                # A point on the cycle
    return max_cycle,point_on_max



def search_quadr(max_heigth):
    max_ab = math.floor(math.exp(max_heigth))      # Recall h(a/b) = log max {|a|,|b|}, so we want to search |a|,|b|<=max_ab 
    max_cycle_length = 0
    for a in range(1,max_ab+1):
        print("Trying a=",a)
        for b in range(1,max_ab+1):         # Suffices to consider b>0, leave sign to be dealt with by a
            if b%1==0:
                print("Trying b=",b)
            if math.gcd(a,b)==1:        # Only compute for fractions in simplest terms, to reduce computation
                f = Henon_quadr_given_ab(a,b)
                longest_cycle,point_on_cycle = get_max_cycle_quadr(f,a,b)
                if longest_cycle>max_cycle_length:
                    print(longest_cycle,"is achieved by a/b =",Fraction(a,b))
                    max_cycle_length = longest_cycle

search_quadr(10)