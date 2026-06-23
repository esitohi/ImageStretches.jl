module ImageStretches

using SpecialFunctions

"""
    ImageStretchFunction{T}

Supertype for all stretch functions.
"""
abstract type ImageStretchFunction{T}
end

export ImageStretchFunction

include("ghs.jl")
export GeneralizedHyperbolicStretch

stretch_image(stretch::ImageStretchFunction, image::AbstractArray) = stretch.(image)

end # module ImageStretches
