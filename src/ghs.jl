"""
    GeneralizedHyperbolicStretch{T<:Real} <: ImageStretchFunction{T}

Implements the generalized hyperbolic stretch function.
Description of the function can be found
[here](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html#__Transformation_equations__).
"""
struct GeneralizedHyperbolicStretch{T<:Real, F<:AbstractFloat} <: ImageStretchFunction{T}
    stretch_factor::T
    D::F
    b::T
    SP::T
    LP::T
    HP::T
    function GeneralizedHyperbolicStretch{T, F}(
        stretch_factor,
        b,
        SP,
        LP,
        HP
    ) where {T, F}
        (0 <= LP < HP) || throw(ArgumentError("LP must be between 0 and HP"))
        (LP < HP <= 1) || throw(ArgumentError("HP must be between LP and 1"))
        (LP <= SP <= HP) || throw(
            ArgumentError("symmetry point must be between lp and hp")
        )
        D = exp(stretch_factor) - 1
        return new(stretch_factor, D, b, SP, LP, HP)
    end
end

function GeneralizedHyperbolicStretch(
    stretch_factor,
    b,
    SP,
    LP = 0,
    HP = 1
)
    args = promote(stretch_factor, b, SP, LP, HP)
    D = exp(stretch_factor) - 1
    return GeneralizedHyperbolicStretch{eltype(args), typeof(D)}(args...)
end

function (ghs::GeneralizedHyperbolicStretch)(x)
    if 0 <= x < ghs.LP
        return NormTi(T1, x, ghs)
    elseif ghs.LP <= x < ghs.SP
        return NormTi(T2, x, ghs)
    elseif ghs.SP <= x < ghs.HP
        return NormTi(T3, x, ghs)
    elseif ghs.HP <= x <= 1
        return NormTi(T4, x, ghs)
    else
        throw(ArgumentError("x must be in [0, 1]"))
    end
end



function BaseLogarithmic(x, ghs::GeneralizedHyperbolicStretch)
    return log(1 + ghs.D*x)
end

function BaseLogarithmicDerivative(x, ghs::GeneralizedHyperbolicStretch)
    return ghs.D/(1 + ghs.D*x)
end


function BaseIntegral(x, ghs::GeneralizedHyperbolicStretch)
    return (1 - (1 - ghs.b*ghs.D*x)^((ghs.b + 1)/ghs.b))/(ghs.D*(ghs.b + 1))
end

function BaseIntegralDerivative(x, ghs::GeneralizedHyperbolicStretch)
    return (1 - ghs.b*ghs.D*x)^(1/ghs.b)
end


function BaseExponential(x, ghs::GeneralizedHyperbolicStretch)
    return 1 - exp(-ghs.D*x)
end

function BaseExponentialDerivative(x, ghs::GeneralizedHyperbolicStretch)
    return ghs.D*exp(-ghs.D*x)
end


function BaseHarmonic(x, ghs::GeneralizedHyperbolicStretch)
    return 1 - (1 + ghs.D*x)^(-1)
end

function BaseHarmonicDerivative(x, ghs::GeneralizedHyperbolicStretch)
    return ghs.D*(1 + ghs.D*x)^(-2)
    
end


function BaseHyperbolic(x, ghs::GeneralizedHyperbolicStretch)
    return 1 - (1 + ghs.b*ghs.D*x)^(-1/ghs.b)
end

function BaseHyperbolicDerivative(x, ghs::GeneralizedHyperbolicStretch)
    return ghs.D*(1 + ghs.b*ghs.D*x)^(-(1 + ghs.b)/ghs.b)
end


function T(x, ghs::GeneralizedHyperbolicStretch)
    if ghs.b == -1
        return BaseLogarithmic(x, ghs)
    elseif ghs.b < 0
        return BaseIntegral(x, ghs)
    elseif ghs.b == 0
        return BaseExponential(x, ghs)
    elseif ghs.b < 1
        return BaseHyperbolic(x, ghs)
    elseif ghs.b == 1
        return BaseHarmonic(x, ghs)
    else
        throw(ArgumentError("b must be in [-1, 1]"))
    end
end

function TDerivative(x, ghs::GeneralizedHyperbolicStretch)

    if ghs.b == -1
        return BaseLogarithmicDerivative(x, ghs)
    elseif ghs.b < 0
        return BaseIntegralDerivative(x, ghs)
    elseif ghs.b == 0
        return BaseExponentialDerivative(x, ghs)
    elseif ghs.b < 1
        return BaseHyperbolicDerivative(x, ghs)
    elseif ghs.b == 1
        return BaseHarmonicDerivative(x, ghs)
    else
        throw(ArgumentError("b must be in [-1, 1]"))
    end
end


function T3(x, ghs::GeneralizedHyperbolicStretch)
    return T(x - ghs.SP, ghs)
end

function T3Derivative(x, ghs::GeneralizedHyperbolicStretch)
    return TDerivative(x - ghs.SP, ghs)
    
end


function T2(x, ghs::GeneralizedHyperbolicStretch)
    return -T(ghs.SP - x, ghs)
end

function T2Derivative(x, ghs::GeneralizedHyperbolicStretch)
    return TDerivative(ghs.SP - x, ghs)
end


function T1(x, ghs::GeneralizedHyperbolicStretch)
    return T2Derivative(x, ghs) * (x - ghs.LP) + T2(ghs.LP, ghs)
end


function T4(x, ghs::GeneralizedHyperbolicStretch)
    return T3Derivative(ghs.HP, ghs) * (x - ghs.HP) + T3(ghs.HP, ghs)
end

function NormTi(Ti::Function, x, ghs::GeneralizedHyperbolicStretch)
    return (Ti(x, ghs) - T1(0, ghs)) / (T4(1, ghs) - T1(0, ghs))
end

