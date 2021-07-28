#
# library of financial models
#
# for now just implementing heston model
#

using QuadGK

#
# underlying asset models (equities)
#

abstract type AbstractModel end

mutable struct HestonModel <: AbstractModel
    S₀::Float64
    r::Float64
    θ::Vector{Float64}
    φ::Function
    
    function HestonModel(S₀::Float64, 
                         r::Float64; 
                         θ=[3.0, 0.1, 0.25, -0.8, 0.08]::Vector{Float64})

        @assert length(θ) == 5 "should be 5 parameters idiot!"

        v₀, v̄, ρ, κ, σ = θ

        ξ(u) = κ - im * σ * ρ * u

        d(u) = sqrt(ξ^2 + σ^2 * (u^2 + im * u))

        A₁(u, t) = (u^2 + im * u) * sinh(d(u) * t / 2)
        A₂(u, t) = d(u) * cosh(d(u) * t / 2) + ξ(u) * sinh(d(u) * t /2)

        A(u, t) = A₁(u, t) / A₂(u, t)

        D(u, t) = log(d(u)) + 0.5 * (κ - d) * t 
                  - log(0.5 * (d(u) + ξ(u) + (d(u) - ξ(u)) * exp(-d(u) * t))) 

        φ(u, t) = exp(im * u * (log(S₀ + r * t) - t * κ * v̄ * ρ / σ) 
                      - (v₀ * A(u, t)) + (2 * κ * v̄ * D(u, t)) / σ^2)

        return new(S₀, r, θ, φ)
    end
    
    HestonModel() = HestonModel(1.0, 0.02)
end

#
# option models
#

abstract type AbstractOption end

mutable struct VanillaCallOption{M <: AbstractModel} <: AbstractOption
    K::Float64
    T::Int
    model::M
end

#
# option pricing functions
# 

function price(option::VanillaCallOption)
    K = option.K
    T = float(option.T) 

    φ = option.model.φ

    I₁ = quadgk(real(exp(-im * u * log(K)) * φ(u - im, T) / (im * u)), 0, Inf)
    I₂ = quadgk(real(exp(-im * u * log(K)) * φ(u, T) / (im * u)), 0, Inf)

    S₀ = option.model.S₀
    r = option.model.r

    return 0.5 * (S₀ - exp(-r * T) * K) + exp(-r * T) / π * (I₁ - K * I₂)
end

