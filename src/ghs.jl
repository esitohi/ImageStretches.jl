#---Preliminary function defintions----------------------------------------------------------------#

function _ghs_T(D, b, x)
    if b == 0
        return -expm1(D*x)                                  # exponential
    elseif b == -1
        return log1p(D*x)                                   # logarithmic
    elseif b == 1
        return 1 - (1 + D*x)^(-1)                           # harmonic
    elseif b < 0
        return (1 - (1 - b*D*x)^((b + 1)/b))/(D*(b + 1))    # integral
    else
        return 1 - (1 + b*D*x)^(-1/b)                       # hyperbolic
    end
end

function _ghs_T_derivative(D, b, x)
    if b == 0
        return D*exp(-D*x)                                  # exponential
    elseif b == -1
        return D/(1 + D*x)                                  # logarithmic
    elseif b == 1
        return D*(1 + D*x)^(-2)                             # harmonic
    elseif b < 0
        return (1 - b*D*x)^(1/b)                            # integral
    else
        return D*(1 + b*D*x)^(-(1 + b)/b)                   # hyperbolic
    end
end

#---GHS type---------------------------------------------------------------------------------------#
"""
    GeneralizedHyperbolicStretch <: ImageStretchFunction

Implements the generalized hyperbolic stretch function.
Description of the function can be found
[here](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html#__Transformation_equations__).
"""
struct GeneralizedHyperbolicStretch <: ImageStretchFunction
    stretch_factor::Float32
    b::Float32
    SP::Float32
    LP::Float32
    HP::Float32
    # Precomputed constants
    _D::Float32     # from stretch factor
    _gT1::Float32   # from LP
    _gT4::Float32   # from HP
    _denom::Float32
    function GeneralizedHyperbolicStretch(stretch_factor, b, SP, LP, HP)
        (0 <= LP < HP) || throw(ArgumentError("LP must be between 0 and HP"))
        (LP < HP <= 1) || throw(ArgumentError("HP must be between LP and 1"))
        (LP <= SP <= HP) || throw(
            ArgumentError("symmetry point must be between LP and HP")
        )
        _D = expm1(stretch_factor)
        _gT1 = -_ghs_T_derivative(_D, b, SP - LP) * -LP - _ghs_T(_D, b, SP - LP)
        _gT4 = _ghs_T_derivative(_D, b, HP - SP) * (1 - HP) - _ghs_T(_D, b, HP - SP)
        _denom = inv(_gT4 - _gT1)
        return new(stretch_factor, b, SP, LP, HP, _D, _gT1, _gT4, _denom)
    end
end

function GeneralizedHyperbolicStretch(stretch_factor, b, SP)
    return GeneralizedHyperbolicStretch(stretch_factor, b, SP, 0, 1)
end

function (ghs::GeneralizedHyperbolicStretch)(x)
    if 0 <= x < ghs.LP
        y = _ghs_T1(ghs, x)
    elseif ghs.LP <= x < ghs.SP
        y = _ghs_T2(ghs, x)
    elseif ghs.SP <= x < ghs.HP
        y = _ghs_T3(ghs, x)
    elseif ghs.HP <= x <= 1
        y = _ghs_T4(ghs, x)
    else
        throw(ArgumentError("x must be in [0, 1]"))
    end
    return (y - ghs._gT1) * ghs._denom
end


base_logarithmic(ghs::GeneralizedHyperbolicStretch, x) = log1p(ghs._D*x)
base_logarithmic_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs._D/(1 + ghs._D*x)


function base_integral(ghs::GeneralizedHyperbolicStretch, x)
    return (1 - (1 - ghs.b*ghs._D*x)^((ghs.b + 1)/ghs.b))/(ghs._D*(ghs.b + 1))
end

function base_integral_derivative(ghs::GeneralizedHyperbolicStretch, x)
    return (1 - ghs.b*ghs._D*x)^(1/ghs.b)
end


base_exponential(ghs::GeneralizedHyperbolicStretch, x) = -expm1(-ghs._D*x)
base_exponential_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs._D*exp(-ghs._D*x)

base_harmonic(ghs::GeneralizedHyperbolicStretch, x) = 1 - (1 + ghs._D*x)^(-1)
base_harmonic_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs._D*(1 + ghs._D*x)^(-2)


function base_hyperbolic(ghs::GeneralizedHyperbolicStretch, x)
    return 1 - (1 + ghs.b*ghs._D*x)^(-1/ghs.b)
end

function base_hyperbolic_derivative(ghs::GeneralizedHyperbolicStretch, x)
    return ghs._D*(1 + ghs.b*ghs._D*x)^(-(1 + ghs.b)/ghs.b)
end


function _ghs_T(ghs::GeneralizedHyperbolicStretch, x)
    if ghs.b == 0
        return base_exponential(ghs, x)
    elseif ghs.b == -1
        return base_logarithmic(ghs, x)
    elseif ghs.b == 1
        return base_harmonic(ghs, x)
    elseif ghs.b < 0
        return base_integral(ghs, x)
    else
        return base_hyperbolic(ghs, x)
    end
end

function _ghs_T_derivative(ghs::GeneralizedHyperbolicStretch, x)
    if ghs.b == -1
        return base_logarithmic_derivative(ghs, x)
    elseif ghs.b < 0
        return base_integral_derivative(ghs, x)
    elseif ghs.b == 0
        return base_exponential_derivative(ghs, x)
    elseif ghs.b == 1
        return base_harmonic_derivative(ghs, x)
    else
        return base_hyperbolic_derivative(ghs, x)
    end
end


_ghs_T3(ghs::GeneralizedHyperbolicStretch, x) = _ghs_T(ghs, x - ghs.SP)
_ghs_T3_derivative(ghs::GeneralizedHyperbolicStretch, x) = _ghs_T_derivative(ghs, x - ghs.SP)

_ghs_T2(ghs::GeneralizedHyperbolicStretch, x) = -_ghs_T(ghs, ghs.SP - x)
_ghs_T2_derivative(ghs::GeneralizedHyperbolicStretch, x) = _ghs_T_derivative(ghs, ghs.SP - x)


function _ghs_T1(ghs::GeneralizedHyperbolicStretch, x)
    return _ghs_T2_derivative(ghs, ghs.LP) * (x - ghs.LP) + _ghs_T2(ghs, ghs.LP)
end


function _ghs_T4(ghs::GeneralizedHyperbolicStretch, x)
    return _ghs_T3_derivative(ghs, ghs.HP) * (x - ghs.HP) + _ghs_T3(ghs, ghs.HP)
end

#---Show methods-----------------------------------------------------------------------------------#

function Base.show(io::IO, ghs::GeneralizedHyperbolicStretch)
    print(io, typeof(ghs), (ghs.stretch_factor, ghs.b, ghs.SP, ghs.LP, ghs.HP))
end

function Base.summary(io::IO, ghs::GeneralizedHyperbolicStretch)
    if ghs.b == 0
        str = "exponential"
    elseif ghs.b == -1
        str = "logarithmic"
    elseif ghs.b == 1
        str = "harmonic"
    elseif ghs.b < 0
        str = "integral"
    else
        str = "hyperbolic"
    end
    print(io, "Generalized hyperbolic stretch (", str, ")")
end

function Base.show(io::IO, ::MIME"text/plain", ghs::GeneralizedHyperbolicStretch)
    summary(io, ghs)
    print(io, ":")
    print(io, '\n', "Stretch factor:       ", ghs.stretch_factor)
    print(io, '\n', "Local intensity:      ", ghs.b)
    print(io, '\n', "Symmetry point:       ", ghs.SP)
    if !iszero(ghs.LP) || !isone(ghs.HP)
        print(io, '\n', "Shadow protection:    ", ghs.LP)
        print(io, '\n', "Highlight protection: ", ghs.HP)
    end
end
