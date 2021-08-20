#
# data generators
#

using LevyProcesses

using Distributions

function volatility_surface(Δs    :: Vector{Tuple{Symbol,Float64}}, 
                            Ts    :: Vector{N},
                            model :: AbstractModel) where {N<:Number}

    surface = Matrix{Float64}(undef, length(Ts), length(Δs))
    
    for (i, T) in enumerate(Ts)
        for (j, (type, K)) in enumerate(Δs)
            if type == :call
                option = VanillaCallOption(K, T, model)
                surface[i,j] = delta(option)
            elseif type == :put
                option = VanillaPutOption(K, T, model)
                surface[i,j] = delta(option)
            else
                error("Unknown option type")
            end
        end
    end

    return surface
end