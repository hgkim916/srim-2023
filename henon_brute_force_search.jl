using Base.Threads
using LinearAlgebra     # needed for norms

function primes_sieve(n)        # Erathostenes - returns all the primes less than or equal to n 
    a = [true for i = 1:n]
    primes = []
    a[1] = false
    for p = 2:n
        if a[p]
            push!(primes,p)
            for i = 2*p:p:n
                a[i] = false
            end
        end
    end
    return primes
end

function Henon_quadr_given_ab(a, b, x_coeff=-1)     # Returns the henon quadratic map
    function f(X)
        x = X[1]
        y = X[2]
        return [y, x_coeff*x + (y^2 + a//b)]
    end
    return f
end

function get_euclidean_bound_quadr(a, b)        # Returns the integer l_inf bound  
    return ceil(Int,1+sqrt(1-a//b))
end

function get_p_adic_bound_quadr(a, b,all_primes)           # Returns the p-adic bound, as the "worst possible denominator" that we could have.
    bound = 1
    for p in all_primes
        if p>b && p > 3
            break
        elseif denominator(b//p)==1      # i.e. p divides b
            # Find exponent of p in fact of b:
            exponent = 0
            while denominator(b//p)==1
                exponent +=1
                b = b//p
            end
            bound = bound * p^(exponent+1)
        elseif p == 2 || p == 3
            bound = bound * p
        end
    end
    return bound
end

function trace_pt(f, eucl_bound, padic_bound, X)    # Gives the orbit of the point X under iteration of f (with escape under the given bounds)
    orbit = Vector{Vector{Rational{Int64}}}()
    while X ∉ orbit
        push!(orbit, X)
        X = f(X)
        x = X[1]
        y = X[2]
        if abs(x) >= eucl_bound || abs(y)>= eucl_bound || denominator(x*padic_bound) > 1 || denominator(y*padic_bound) > 1
            # We iterate outside to somewhere we can guarantee we're not preperiodic
            return nothing
        end
    end
    return orbit
end

function get_max_cycle_quadr(f, a, b,all_primes)    # return the longest cycle of the quadratic Henon map f given by a and b, as well as a point on the cycle
    eucl_bound = get_euclidean_bound_quadr(a, b)
    padic_bound = get_p_adic_bound_quadr(a, b, all_primes)
    max_cycle = 0
    point_on_max = []
    checked_points = []
    
    # search all (x1/x2,y1/y2) within the bounds by multiplying by the "worst possible denominator" and searching in the integers
    search_space = eucl_bound * padic_bound     
    for x in -search_space:search_space
        for y in -search_space:search_space
            X= [x//padic_bound, y//padic_bound]
            if X ∉ checked_points       # only check once, to shorten computation
                orbit = trace_pt(f, eucl_bound, padic_bound, X)
                if orbit ≠ nothing
                    append!(checked_points, orbit)
                    if length(orbit) > max_cycle    # longest cycle so far, save
                        max_cycle = length(orbit)
                        point_on_max = X
                    end
                end
            end
        end
    end
    return max_cycle, point_on_max
end

function search_quadr(max_height)           # searches among all the quadratics desired with height at most max_height
    max_ab = floor(Int,exp(max_height))     # Recall h(a,b) = log max(|a|,|b|), we invert this to get a range for |a| and |b|
    println("height ",max_height," means search to modulus ",max_ab)
    max_cycle_length = Threads.Atomic{Int}(0)
    all_primes = primes_sieve(max_ab)
    for a in 0:-1:-max_ab     
        Threads.@threads for b in 1:max_ab
            if gcd(a, b) == 1
                f = Henon_quadr_given_ab(a, b)
                longest_cycle, point_on_cycle = get_max_cycle_quadr(f, a, b,all_primes)
                Threads.atomic_max!(max_cycle_length, longest_cycle)
                if longest_cycle >= max_cycle_length[]       # better than anything so far
                    println(longest_cycle, " is achieved by a/b = ", a//b," starting at [",point_on_cycle[1],",",point_on_cycle[2],"]")
                    orbit = trace_pt(f,get_euclidean_bound_quadr(a,b),get_p_adic_bound_quadr(a,b,all_primes),point_on_cycle)
                    println("Orbit achieving this is: ",orbit)
                end
            end
        end
        println("done with a=",a)
    end
    println("done!")
end

search_quadr(5)