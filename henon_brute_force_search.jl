using Base.Threads
using LinearAlgebra

function Henon_quadr_given_ab(a, b, x_coeff=-1)
    function f(X)
        x = X[1]
        y = X[2]
        return [y, x_coeff*x + (y^2 + a//b)]
    end
    return f
end

function get_euclidean_bound_quadr(a, b)
    return 2*(1+max(1,abs(a//b)))
end

function get_p_adic_bound_quadr(a, b)
    return b
end

function trace_pt(f, eucl_bound, padic_bound, X)
    orbit = Vector{Vector{Rational{Int64}}}()
    while X ∉ orbit
        push!(orbit, X)
        X = f(X)
        x = X[1]
        y = X[2]
        if norm(X) > eucl_bound || denominator(x*padic_bound) > 1 || denominator(y*padic_bound) > 1
            return nothing
        end
    end
    return orbit
end

function get_max_cycle_quadr(f, a, b)
    eucl_bound = get_euclidean_bound_quadr(a, b)
    padic_bound = get_p_adic_bound_quadr(a, b)
    max_cycle = 0
    point_on_max = []
    checked_points = []
    
    search_space = eucl_bound * padic_bound
    for x in -search_space:search_space
        for y in -search_space:search_space
            X= [x//padic_bound, y//padic_bound]
            if X ∉ checked_points
                orbit = trace_pt(f, eucl_bound, padic_bound, X)
                if orbit ≠ nothing
                    append!(checked_points, orbit)
                    if length(orbit) > max_cycle
                        max_cycle = length(orbit)
                        point_on_max = X
                    end
                end
            end
        end
    end
    return max_cycle, point_on_max
end

function search_quadr(max_height)
    max_ab = floor(Int,exp(max_height))
    println("height ",max_height," means search to modulus ",max_ab)
    max_cycle_length = Threads.Atomic{Int}(0)
    for a in 1:max_ab
        Threads.@threads for b in 1:max_ab
            if gcd(a, b) == 1
                f = Henon_quadr_given_ab(a, b)
                longest_cycle, point_on_cycle = get_max_cycle_quadr(f, a, b)
                if longest_cycle > max_cycle_length[]
                    println(longest_cycle, " is achieved by a/b = ", a//b," starting at [",point_on_cycle[1],",",point_on_cycle[2],"]")
                    orbit = trace_pt(f,get_euclidean_bound_quadr(a,b),get_p_adic_bound_quadr(a,b),point_on_cycle)
                    println("Orbit achieving this is: ",orbit)
                    max_cycle_length[] = longest_cycle
                end
            end
        end
        println("done with a=",a)
    end
    println("done!")
end

search_quadr(5)