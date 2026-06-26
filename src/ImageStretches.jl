module ImageStretches

using SpecialFunctions

"""
    ImageStretchFunction

Supertype for all stretch functions.

Instances of stretch functions should be functors (in other words, they can be called as
functions) which can be applied to a real argument.
The functor should try to return the same type as the input.
It is recommended that all field types be `Float32`, since this is the highest precision
commonly used in image processing.
"""
abstract type ImageStretchFunction
end

export ImageStretchFunction

include("ghs.jl")
export GeneralizedHyperbolicStretch

function stretch_image(stretch::ImageStretchFunction, image::AbstractArray)
    result = similar(image)
    # @inbounds makes this slightly faster than map(), broadcasting, or a comprehension
    @inbounds for i in eachindex(image)
        result[i] = stretch(image[i])
    end
    return result
end

end # module ImageStretches
