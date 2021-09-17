module RationalApproximations

using LsqFit

export value, fit, RationalPoly

struct RationalPoly
    numerator::Vector{Float64}
    denominator::Vector{Float64}
end

to_num_den(p::Vector{Float64}, deg_num::Int) = (p[1:deg_num+1], p[deg_num+2:end])

function RationalPoly(p::Vector{Float64}, degree_numerator::Int)
    num, den = to_num_den(p, degree_numerator)
    return RationalPoly(num, den)
end

function eval_poly(x, coef::Vector{Float64})
    out = coef[1] .+ zero(x)

    for (k, c) in enumerate(coef[2:end])
        out += c*x.^k
    end

    return out
end

function value(x, p::Vector{Float64}, deg_num::Int) 
    num, den = to_num_den(p, deg_num)
    return eval_poly(x, num) ./ eval_poly(x, den)
end
value(x, p::RationalPoly) = eval_poly(x, p.numerator) ./ eval_poly(x, p.denominator)


# computes the gradient with respect to the coefficients
function grad(x, p::Vector{Float64}, deg_num::Int)
    num, den = to_num_den(p, deg_num)

    x_all_num = [1; [x^(k-1) for k in 2:length(num)]]
    grad_num = x_all_num / eval_poly(x, den)

    x_all_den = [1; [x^(k-1) for k in 2:length(den)]]
    grad_den = -x_all_den * eval_poly(x, num) / eval_poly(x, den)^2

    return [grad_num; grad_den]
end

# A lot of allocations, but w/e, should be fine for now
function jac(t, p, deg_num::Int)
    J = zeros(length(t), length(p))

    for (i, v) in enumerate(t)
        J[i,:] = grad(v, p, deg_num)
    end

    return J
end

function fit(x, y, degree_numerator::Int, degree_denominator::Int)
    @assert degree_numerator >= 0 && degree_denominator >= 0
    coef = ones(degree_numerator+degree_denominator+2)

    m(t, p) = value(t, p, degree_numerator)
    j(t, p) = jac(t, p, degree_numerator)

    fit = curve_fit(m, j, x, y, coef)

    return RationalPoly(fit.param, degree_numerator)
end

end
