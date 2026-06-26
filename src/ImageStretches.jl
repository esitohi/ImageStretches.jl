module ImageStretches

using SpecialFunctions

"""
    ImageStretchFunction{T}

Supertype for all stretch functions.

Instances of stretch functions should be functors (in other words, they can be called as
functions) which can be applied to a real argument.
The functor should avoid changing the type of the input if possible.
"""
abstract type ImageStretchFunction{T}
end

export ImageStretchFunction

include("ghs.jl")
export GeneralizedHyperbolicStretch

stretch_image(stretch::ImageStretchFunction, image::AbstractArray) = stretch.(image)

end # module ImageStretches
