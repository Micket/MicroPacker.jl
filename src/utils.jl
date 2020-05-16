import Distributions: Uniform, Weibull

uniform_dist(l, u) = return () -> rand(Uniform(l, u))
weibull_dist(a, mu) = return () -> rand(Weibull(a, mu))

function parse_dist(arg)
    # Parses input string for given distribution.
    # Returns a distribution, and the average
    d, params = split(arg, ':')
    params = [parse(Float64, p) for p in split(params, ',')]
    if d == "U"
        return uniform_dist(params...), sum(params)/2
    elseif d == "W"
        a, mu = params
        return weibull_dist(a, mu), mu
    elseif d == "C"
        return () -> params[1], params[1]
    else
        throw(DomainError(d, "Unrecognized distribution function"))
    end
end

mass_fraction(vol_frac, ρ_i, ρ_m) = vol_frac*ρ_i / (vol_frac*ρ_i + (1-vol_frac)*ρ_m)
