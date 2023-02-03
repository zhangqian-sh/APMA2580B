function draw_characteristics_of_burgers_equation_with_sin_initial(N::Int)
    """
    uₜ + (u^2 / 2)ₓ = 0
    u(0, x) = sin(x)
    x ∈ [-π, π]
    """
    u₀(x) = sin(x)
    du₀(x) = cos(x)
    list_x₀ = LinRange(0, π, N + 1)
    list_k = u₀.(list_x₀)

    function evaluate(t::Float64, x::Float64, tol::Float64=1e-16, maxiter::Int=100)
        function evaluate_in_0_to_pi(t::Float64, x::Float64, tol::Float64=1e-16, maxiter::Int=100)
            # find nearest characteristics
            d = zeros(Float64, N + 1)
            for i = 1:N+1
                xₚ = list_k[i] * t + list_x₀[i]
                d[i] = abs(xₚ - x)
            end
            nearest_characteristics_idx = argmin(d)

            # use Newton iteration with the nearest characteristics as the initial guess
            xₜ = list_x₀[nearest_characteristics_idx]
            F(x₀) = x₀ + t * u₀(x₀) - x
            dF(x₀) = 1 + t * du₀(x₀)
            for _iter = 1:maxiter
                xₜ = xₜ - F(xₜ) / dF(xₜ)
                if abs(F(xₜ)) < tol
                    break
                end
            end

            # use the characteristics through x to evaluate u
            u = u₀(xₜ)
        end

        # put x into [-π, π]
        x -= floor((x + π) / 2π) * 2π
        if x > 0
            evaluate_in_0_to_pi(t, x, tol, maxiter)
        else
            -evaluate_in_0_to_pi(t, -x, tol, maxiter)
        end
    end
end

α, β = 0.5, 1
t, x = 2.0, 4.1
N_characteristics = 200
evaluate = draw_characteristics_of_burgers_equation_with_sin_initial(N_characteristics)
v = evaluate(β * t, x - α * t)
u = β * v + α
println(u)
