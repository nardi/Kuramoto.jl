module Kuramoto

export Kuramoto, kuramoto

using DynamicalSystemsBase, Distributions, Random

"""
    Kuramoto.system(ω::AbstractArray;
                    θ0=rand(Uniform(-π, π), size(ω)), K=0)
    Kuramoto.system(N::Int, ω_dist::Sampleable;
                θ0_dist=Uniform(-π, π), θ0=rand(θ0_dist, N), K=0)

    Alias: kuramoto

    Constructs an instance of the Kuramoto system. Constructs N entities,
    each of which moves at one of the natural velocities ω, and `attracts'
    the others with a `force' which is proportional to the constant K and
    decreases with (circle) distance. θ0 gives the starting positions for
    each of the entities. 
"""
function system(ω::AbstractArray;
                θ0=rand(Uniform(-π, π), size(ω)), K=0)
    N = size(ω)[1]
    function kuramoto_eom!(dθ, θ, (ω, K), t)
        C = K / N
        c = @. exp(θ*im)
        S = sum(c)
        dθ .= @. ω + C*imag(S * conj(c))
    end
    
    ContinuousDynamicalSystem(kuramoto_eom!, θ0, (ω, K))    
end

function system(N::Int, ω_dist::Sampleable;
                θ0_dist=Uniform(-π, π), θ0=rand(θ0_dist, N), K=0)
    ω = rand(ω_dist, N);
    system(ω; θ0=θ0, K=K)
end

kuramoto = system

"""
    Kuramoto.critical(ω_dist = Normal(0, 1))

    Calculates the theoretical critical value of K for
    a given distribution of natural velocities.
"""
critical(ω_dist = Normal(0, 1)) = 2/(π*pdf(ω_dist, 0))

"""
    Kuramoto.order(trj)

    Calculates the order parameter (the concentration
    of the angles) for a trajectory (should be a tall
    array or Dataset).
"""
function order(trj)
    N = size(trj)[2]
    map(θ -> abs(sum(exp.(im*θ)/N)), trj)
end

end # module
