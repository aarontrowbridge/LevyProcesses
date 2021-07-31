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
    S :: Float64
    r :: Float64
    θ :: Vector{Float64}
    φ :: Function
    
    function HestonModel(S :: Float64, 
                         r :: Float64; 
                         θ=[3.0, 0.1, 0.25, -0.8, 0.08]::Vector{Float64})

        @assert length(θ) == 5 "should be 5 parameters idiot!"

        κ, v̄, σ, ρ, v₀ = θ

        ξ(u) = κ - im * σ * ρ * u

        d(u) = sqrt(ξ(u)^2 + σ^2 * (u^2 + im * u))

        A₁(u, t) = (u^2 + im * u) * sinh(d(u) * t / 2)
        A₂(u, t) = d(u) * cosh(d(u) * t / 2) + ξ(u) * sinh(d(u) * t / 2)

        A(u, t) = A₁(u, t) / A₂(u, t)

        D(u, t) = log(d(u)) + 0.5 * (κ - d(u)) * t 
                  - log(0.5 * (d(u) + ξ(u) + (d(u) - ξ(u)) * exp(-d(u) * t))) 

        φ(u, t) = exp(im * u * (log(S) + r * t - t * κ * v̄ * ρ / σ) 
                      - v₀ * A(u, t) + 2 * κ * v̄ / σ^2 * D(u, t))

        return new(S, r, θ, φ)
    end
    
    HestonModel() = HestonModel(1.0, 0.02)
end

#
# option models
#

abstract type AbstractOption end

mutable struct VanillaCallOption{M <: AbstractModel} <: AbstractOption
    K     :: Float64
    T     :: Number
    model :: M
end

mutable struct VanillaPutOption{M <: AbstractModel} <: AbstractOption
    K     :: Float64
    T     :: Number
    model :: M 
end

#
# option pricing functions
# 

function price(option::VanillaCallOption)
    K = option.K
    T = float(option.T) 

    φ = option.model.φ

    f₁(u) = real(exp(-im * u * log(K)) / (im * u) * φ(u - im, T))
    f₂(u) = real(exp(-im * u * log(K)) / (im * u) * φ(u, T))

    I₁, err1 = quadgk(f₁, 0, 10)
    I₂, err2 = quadgk(f₂, 0, 10)

    S = option.model.S
    r = option.model.r

    return 0.5 * (S - exp(-r * T) * K) + exp(-r * T) / π * (I₁ - K * I₂)
end

function price(option::VanillaPutOption)
    call = VanillaCallOption(option.K, option.T, option.model)
    return price(call) - option.model.S + option.K * exp(-option.model.r * option.T)
end

function delta(option::VanillaCallOption)
    K = option.K
    T = float(option.T)

    φ = option.model.φ

    f₁(u) = real(exp(-im * u * log(K)) * (u - im) / u * φ(u - im, T))
    f₂(u) = real(exp(-im * u * log(K)) * φ(u, T))

    ∂I₁, err1 = quadgk(f₁, 0, 10)
    ∂I₂, err2 = quadgk(f₂, 0, 10)

    S = option.model.S
    r = option.model.r

    return 0.5 + exp(-r * T) / (π * S) * (∂I₁ - K * ∂I₂)
end

function delta(option::VanillaPutOption)
    call = VanillaCallOption(option.K, option.T, option.model)
    return delta(call) - 1
end
