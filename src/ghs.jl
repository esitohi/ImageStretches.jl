"""
    GeneralizedHyperbolicStretch{T<:Real}

Implements the generalized hyperbolic stretch function.
Description of the function can be found
[here](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html#__Transformation_equations__).
"""
struct GeneralizedHyperbolicStretch{T<:Real} <: ImageStretchFunction{T}
    stretch_factor::T
    local_intensity::T
    symmetry_point::T
    lp::T
    hp::T
    function GeneralizedHyperbolicStretch{T}(
        stretch_factor,
        local_intensity,
        symmetry_point,
        lp,
        hp
    ) where T
        (0 <= lp < hp) || throw(ArgumentError("lp must be between 0 and hp"))
        (lp < hp <= 1) || throw(ArgumentError("hp must be between lp and 1"))
        (lp <= symmetry_point <= hp) || throw(
            ArgumentError("symmetry point must be between lp and hp")
        )
        return new(stretch_factor, local_intensity, symmetry_point, lp, hp)
    end
end

function GeneralizedHyperbolicStretch(
    stretch_factor,
    local_intensity,
    symmetry_point,
    lp = false,
    hp = true)
    args = promote(stretch_factor, local_intensity, symmetry_point, lp, hp)
    return GeneralizedHyperbolicStretch{eltype(args)}(args...)
end

#=
function (ghs::GeneralizedHyperbolicStretch)(x)
    
end
=#
