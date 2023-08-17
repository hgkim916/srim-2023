import math

#define the polys that we need:
def f2(x):
    return (22-9*x+x**2)/2

def f3(x):
    return (-66+89*x-18*x**2+x**3)/6

def f4_1(x):
    return (312-374*x+155*x**2-22*x**3+x**4)/24

def f4_2(x):
    return (552-506*x+167*x**2-22*x**3+x**4)/24

def f4_3(x):
    return (408-440*x+161*x**2-22*x**3+x**4)/12

def f5(x):
    return (-3600+5794*x-2485*x**2+445*x**3-35*x**4+x**5)/120

def f6(x):
    return (30960-45060*x+25504*x**2-6375*x**3+775*x**4-45*x**5+x**6)/720

def f7(x):
    return (-151200+373764*x-258104*x**2+83629*x**3-14000*x**4+1246*x**5-56*x**6+x**7)/5040

def escape(x,f,escape_radius): # checks if x escapes the escape radius under iteration of f
    x_init = x
    iter = 500000       # adjust: higher takes longer, but means more certain that it's indeed preperiodic
    for i in range(iter):       
        x = f(x)
        if x_init.is_integer():         # this just adds precision if we're working with integers, in the hope or reducing rounding errors
            x=int(x)
        if abs(x) > escape_radius:          # the iterate has left the escape radius, and hence escapes to infinity, so can stop the check
            # print("escapes to infinty, operation aborted")
            return True
    # print(x_init, "reaches an end value of",x)
    return False        # the iterates have not yet escaped after the large number of iterations, so we conclude likely preperiodic

def check_rat_preper(eucl_bd,pbound,f):     # looks for rational preperiodic points
                                            # eucl_bound is the maximum euclidean norm such a point can have
                                            # pbound is the maximum denominator a preperiodic point can have (obtained from the p-adic norm bounds)
    minint = - eucl_bd
    maxint =  eucl_bd
    print("The rational preperiodic points are:")
    for x in range(minint*pbound,maxint*pbound+1,1):      # x is the numerator, pbound is the denominator, and we iterate from minint to maxint in step 1/pbound
        # escape(x/pbound,f)
        if not escape(x/pbound,f,eucl_bd):
            if int(x/pbound)*pbound == int(x):          # print integers as integers
                print(int(x/pbound))
            else:                                       # print rationals as fractions
                print(int(x/math.gcd(x,pbound)),"/",int(pbound/math.gcd(x,pbound)))

def automated(coeff,denom):    # Automates more bits of the check, only receiving as inputs the coefficients of the numerator poly and the (gcd) denominator
    d = len(coeff)-1
    if coeff[d] != 1:
        p_bound = int(input("Could not determine p-adic bound for the polynomial. Please enter: "))
    else:
        p_bound = 1
        
    max_coeff = max([abs(coeff[i]) for i in range(d)])
    Eucl_bound = math.ceil(max(1,abs((1/coeff[d])*(abs(denom)+d*max_coeff))))

    def f(x):
        result = 0
        for i in range(d+1):
            result+=coeff[i]*(x**i)
        # print("x =", x," gives f(x) =",result)
        return result/denom
    


    check_rat_preper(Eucl_bound,p_bound,f)

# CHECK FOR RATIONAL PREPERIODIC POINTS - MORE AUTOMATED


# automated([22,-9,1],2)

# automated([14,-7,1],2)

# automated([10,-6,1],1)

# automated([8,-5,1],2)

# automated([20,-15,3],2)

# automated([4,-3,1],2)

# automated([7,-5,1],1)

# automated([16,-13,3],2)

# automated([12,-11,3],2)

# automated([5,-4,1],1)