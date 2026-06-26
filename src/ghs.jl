"""
    GeneralizedHyperbolicStretch <: ImageStretchFunction{Float32}

Implements the generalized hyperbolic stretch function.
Description of the function can be found
[here](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html#__Transformation_equations__).
"""
struct GeneralizedHyperbolicStretch <: ImageStretchFunction{Float32}
    stretch_factor::Float32
    D::Float32    # precomputed
    b::Float32
    SP::Float32
    LP::Float32
    HP::Float32
    function GeneralizedHyperbolicStretch(stretch_factor, b, SP, LP, HP)
        (0 <= LP < HP) || throw(ArgumentError("LP must be between 0 and HP"))
        (LP < HP <= 1) || throw(ArgumentError("HP must be between LP and 1"))
        (LP <= SP <= HP) || throw(
            ArgumentError("symmetry point must be between LP and HP")
        )
        return new(stretch_factor, expm1(stretch_factor), b, SP, LP, HP)
    end
end

function GeneralizedHyperbolicStretch(stretch_factor, b, SP)
    return GeneralizedHyperbolicStretch(stretch_factor, b, SP, 0, 1)
end

function (ghs::GeneralizedHyperbolicStretch)(x)
    if 0 <= x < ghs.LP
        return _ghs_NormTi(_ghs_T1, ghs, x)
    elseif ghs.LP <= x < ghs.SP
        return _ghs_NormTi(_ghs_T2, ghs, x)
    elseif ghs.SP <= x < ghs.HP
        return _ghs_NormTi(_ghs_T3, ghs, x)
    elseif ghs.HP <= x <= 1
        return _ghs_NormTi(_ghs_T4, ghs, x)
    else
        throw(ArgumentError("x must be in [0, 1]"))
    end
end


base_logarithmic(ghs::GeneralizedHyperbolicStretch, x) = log1p(ghs.D*x)
base_logarithmic_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs.D/(1 + ghs.D*x)


function base_integral(ghs::GeneralizedHyperbolicStretch, x)
    return (1 - (1 - ghs.b*ghs.D*x)^((ghs.b + 1)/ghs.b))/(ghs.D*(ghs.b + 1))
end

function base_integral_derivative(ghs::GeneralizedHyperbolicStretch, x)
    return (1 - ghs.b*ghs.D*x)^(1/ghs.b)
end


base_exponential(ghs::GeneralizedHyperbolicStretch, x) = -expm1(-ghs.D*x)
base_exponential_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs.D*exp(-ghs.D*x)

base_harmonic(ghs::GeneralizedHyperbolicStretch, x) = 1 - (1 + ghs.D*x)^(-1)
base_harmonic_derivative(ghs::GeneralizedHyperbolicStretch, x) = ghs.D*(1 + ghs.D*x)^(-2)


function base_hyperbolic(ghs::GeneralizedHyperbolicStretch, x)
    return 1 - (1 + ghs.b*ghs.D*x)^(-1/ghs.b)
end

function base_hyperbolic_derivative(ghs::GeneralizedHyperbolicStretch, x)
    return ghs.D*(1 + ghs.b*ghs.D*x)^(-(1 + ghs.b)/ghs.b)
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
    return _ghs_T2_derivative(ghs, x) * (x - ghs.LP) + _ghs_T2(ghs, ghs.LP)
end


function _ghs_T4(ghs::GeneralizedHyperbolicStretch, x)
    return _ghs_T3_derivative(ghs, ghs.HP) * (x - ghs.HP) + _ghs_T3(ghs, ghs.HP)
end

function _ghs_NormTi(Ti::Function, ghs::GeneralizedHyperbolicStretch, x)
    return (Ti(ghs, x) - _ghs_T1(ghs, 0)) / (_ghs_T4(ghs, 1) - _ghs_T1(ghs, 0))
end

