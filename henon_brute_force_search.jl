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

# QUADRATIC STUFF

function Henon_quadr_given_ab(a, b, x_coeff=-1)     # Returns the henon quadratic map
    function f(X)
        x = X[1]
        y = X[2]
        return [y, x_coeff*x + (y^2 + a//b)]
    end
    return f
end

function get_euclidean_bound_quadr(a, b)        # Returns the integer l_inf bound  
    if a//b>1       # in this case, no periodic points can exist
        return 0
    else
        return ceil(Int,1+sqrt(1-a//b))
    end
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
    Threads.@threads for a in 0:-1:-max_ab
        Threads.@threads for b in 1:max_ab
            if gcd(a, b) == 1
                f = Henon_quadr_given_ab(a, b)
                longest_cycle, point_on_cycle = get_max_cycle_quadr(f, a, b,all_primes)
                Threads.atomic_max!(max_cycle_length, longest_cycle)
                if longest_cycle >= max_cycle_length[]       # better than anything so far
                    orbit = trace_pt(f,get_euclidean_bound_quadr(a,b),get_p_adic_bound_quadr(a,b,all_primes),point_on_cycle)
                    if longest_cycle>= max_cycle_length[]       # check again, so that it hasn't been updated because of threading
                        println(longest_cycle, " is achieved by a/b = ", a//b," starting at [",point_on_cycle[1],",",point_on_cycle[2],"]\n Orbit achieving this is: ",orbit,"\n Maximum now is ",max_cycle_length[])
                    end
                end
                println("done with a/b=",a//b)
                if a!=0
                    # We will run the -a case from here, so that the order in which we look at a is 0,1,-1,2,-2,... - I think this makes more sense.
                    a = -a
                    f = Henon_quadr_given_ab(a, b)
                    longest_cycle, point_on_cycle = get_max_cycle_quadr(f, a, b,all_primes)
                    Threads.atomic_max!(max_cycle_length, longest_cycle)
                    if longest_cycle >= max_cycle_length[]       # better than anything so far
                        orbit = trace_pt(f,get_euclidean_bound_quadr(a,b),get_p_adic_bound_quadr(a,b,all_primes),point_on_cycle)
                        if longest_cycle>= max_cycle_length[]       # check again, so that it hasn't been updated because of threading
                            println(longest_cycle, " is achieved by a/b = ", a//b," starting at [",point_on_cycle[1],",",point_on_cycle[2],"]\n Orbit achieving this is: ",orbit,"\n Maximum now is ",max_cycle_length[])
                        end
                    end
                    println("done with a/b=",a//b)
                    a=-a
                end
            end
        end
    end
    println("done!")
end

# GENERAL POLYNOMIALS

function Henon_general_poly(as,bs,x_coeff = -1)         # Returns the Henon map (x,y)->(y,x_coeff*x+sum[(as[i]/bs[i])*y^i])
    function f(X)
        x = X[1]
        y = X[2]
        p = 0
        for i in eachindex(as)
            p += (as[i]//bs[i])*(y^i)
        end
        return [y, x_coeff*x+p]
    end
    return f    
end

function get_euclidean_bound_general(as, bs)        # Returns the integer l_inf bound  
    A = maximum(abs(as[i]//bs[i]) for i in eachindex(as[1:1:length(as)]))
    R = (2+A)//(abs(as[length(as)]//bs[length(bs)]))
    return max(1,ceil(R))
end

function get_p_adic_bound_general(as, bs,all_primes)           # Returns the p-adic bound, as the "worst possible denominator" that we could have.
    bound = 1
    # for p in all_primes
    #     if p>b && p > 3
    #         break
    #     elseif denominator(b//p)==1      # i.e. p divides b
    #         # Find exponent of p in fact of b:
    #         exponent = 0
    #         while denominator(b//p)==1
    #             exponent +=1
    #             b = b//p
    #         end
    #         bound = bound * p^(exponent+1)
    #     elseif p == 2 || p == 3
    #         bound = bound * p
    #     end
    # end
    return bound
end

function get_max_cycle_general(f, as, bs,all_primes)    # return the longest cycle of the quadratic Henon map f given by as and bs
    eucl_bound = get_euclidean_bound_general(as, bs)
    padic_bound = get_p_adic_bound_general(as, bs, all_primes)
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

# function print_poly(as,bs)
#     d = length(as)
#     string = ""
#     for i in eachindex(as)
#         if as[i] != 0
#             if i == 1
#                 if as[i]//bs[i] == 1
#                     string +="+x"
#                 elseif poly_v[i] == -1
#                     string +="-x"
#                 else
#                     string+="x"
#                 end
#             elseif i == 2
#                 if poly_v[i] == 1
#                     print("+x^$i", end="")
#                 elseif poly_v[i] == -1
#                     print("-x^$i", end="")
#                 else
#                     print("$(poly_v[i])x^$i", end="")
#             else
#                 if poly_v[i] == 1
#                     print("+x^$(i)", end="")
#                 elseif poly_v[i] == -1
#                     print("-x^$(i)", end="")
#                 else
#                     print("$(poly_v[i])x^$(i)", end="")
#             end
#         end
#     end
# end

function search_general(max_height,d)           # searches among all the polys of degree d with height at most max_height
    max_ab = floor(Int,exp(max_height))     # Recall h(a,b) = log max(|a|,|b|), we invert this to get a range for |a| and |b|
    println("height ",max_height," means search to modulus ",max_ab)
    max_cycle_length = Threads.Atomic{Int}(0)
    all_primes = primes_sieve(max_ab)

    search_space_a = [i for i in Iterators.product(ntuple(_ -> -max_ab:max_ab,d+1)...)]
    search_space_b = [i for i in Iterators.product(ntuple(_ -> 0:max_ab,d+1)...)]
    Threads.@threads for as in search_space_a
        Threads.@threads for bs in search_space_b
            if max([gcd(as[i],bs[i]) for i in eachindex(as)])==1
                f = Henon_general_poly(as, bs)
                longest_cycle, point_on_cycle = get_max_cycle_general(f,as,bs,all_primes)
                Threads.atomic_max!(max_cycle_length, longest_cycle)
                if longest_cycle >= max_cycle_length[]       # better than anything so far
                    orbit = trace_pt(f,get_euclidean_bound_general(as,bs),get_p_adic_bound_general(as,bs,all_primes),point_on_cycle)
                    if longest_cycle>= max_cycle_length[]       # check again, so that it hasn't been updated because of threading
                        println(longest_cycle, " is achieved by as=",as," and bs=",bs,"\n    Orbit achieving this is: ",orbit,"\n    Maximum now is ",max_cycle_length[])
                    end
                end
            end
        end
    end
    println("done!")
end


# search_quadr(5)
search_general(2,3)