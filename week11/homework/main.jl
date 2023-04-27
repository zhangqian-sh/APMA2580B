using DelimitedFiles
using QuadGK
using Statistics
using StatProfilerHTML

include("reconstruction.jl")
include("burgers_exact_solver.jl")

# Numerical Flux
# Here f(u) = u^2 / 2, f'(u) = u.
f(u)::Float64 = u^2 / 2
g(u)::Float64 = u^2 / 2
# debug: check linear equation
# f(u) = u

# const scheme = ARGS[1]
t_end = parse(Float64, ARGS[1])
const λ_max = 0.6 # max CFL number

# flux splitting
function LaxFriedrichsFlux(uₗ, uᵣ, α)
    (f(uₗ) + f(uᵣ) - α * (uᵣ - uₗ)) / 2
end

function f⁺(uₗ, α)
    (f(uₗ) + α * uₗ) / 2
end

function f⁻(uᵣ, α)
    (f(uᵣ) - α * uᵣ) / 2
end

function g⁺(uₗ, α)
    (g(uₗ) + α * uₗ) / 2
end

function g⁻(uᵣ, α)
    (g(uᵣ) - α * uᵣ) / 2
end

flux = LaxFriedrichsFlux

function solve(u)
    quit_loop = false
    t = 0
    N_y, N_x = size(u)
    Δx = 4π / N_x
    Δy = 4π / N_y
    Δt = λ_max * Δx^(5 / 3)
    λ_x = Δt / Δx
    λ_y = Δt / Δy
    stencil_size = 3

    fₚ, fₙ = zeros(Float64, (N_x)), zeros(Float64, (N_x))
    f̂⁺, f̂⁻ = zeros(Float64, (N_y, N_x + 1)), zeros(Float64, (N_y, N_x + 1))
    f̂ = zeros(Float64, (N_y, N_x + 1))
    gₚ, gₙ = zeros(Float64, (N_y)), zeros(Float64, (N_y))
    ĝ⁺, ĝ⁻ = zeros(Float64, (N_y + 1, N_x)), zeros(Float64, (N_y + 1, N_x))
    ĝ = zeros(Float64, (N_y + 1, N_x))

    v = zeros(Float64, (N_y, N_x))

    periodic_x(j) = periodic_bc(j, N_x)
    periodic_y(i) = periodic_bc(i, N_y)
    reconstruct = WENO5_reconstruct

    # single forward step
    function Euler_forward(u)::Array{Float64}
        # reconstruct flux splitting of f̂
        for i = 1:N_y
            α = maximum(abs.(u[i, :]))
            for j = 1:N_x
                fₚ[j], fₙ[j] = f⁺(u[i, j], α), f⁻(u[i, j], α)
            end
            for j = 1:N_x
                f̂⁻[i, j+1] = reconstruct(
                    fₚ[periodic_x(j - 2)], fₚ[periodic_x(j - 1)], fₚ[periodic_x(j)], fₚ[periodic_x(j + 1)], fₚ[periodic_x(j + 2)]
                )
                f̂⁺[i, j] = reconstruct(
                    fₙ[periodic_x(j + 2)], fₙ[periodic_x(j + 1)], fₙ[periodic_x(j)], fₙ[periodic_x(j - 1)], fₙ[periodic_x(j - 2)]
                )
            end
            f̂⁻[i, 1] = f̂⁻[i, end]
            f̂⁺[i, end] = f̂⁺[i, 1]
            # monotone flux: Lax-Friedrichs
            for j = 1:N_x+1
                f̂[i, j] = f̂⁺[i, j] + f̂⁻[i, j]
            end
        end

        # reconstruct flux splitting of ĝ
        for j = 1:N_x
            α = maximum(abs.(u[:, j]))
            for i = 1:N_y
                gₚ[i], gₙ[i] = g⁺(u[i, j], α), g⁻(u[i, j], α)
            end
            for i = 1:N_y
                ĝ⁻[i+1, j] = reconstruct(
                    gₚ[periodic_y(i - 2)], gₚ[periodic_y(i - 1)], gₚ[periodic_y(i)], gₚ[periodic_y(i + 1)], gₚ[periodic_y(i + 2)]
                )
                ĝ⁺[i, j] = reconstruct(
                    gₙ[periodic_y(i + 2)], gₙ[periodic_y(i + 1)], gₙ[periodic_y(i)], gₙ[periodic_y(i - 1)], gₙ[periodic_y(i - 2)]
                )
            end
            ĝ⁻[1, j] = ĝ⁻[end, j]
            ĝ⁺[end, j] = ĝ⁺[1, j]
            # monotone flux: Lax-Friedrichs
            for i = 1:N_x+1
                ĝ[i, j] = ĝ⁺[i, j] + ĝ⁻[i, j]
            end
        end

        # Euler forward step
        v = u - λ_x * (f̂[:, 2:end] - f̂[:, 1:end-1]) - λ_y * (ĝ[2:end, :] - ĝ[1:end-1, :])
    end

    function SSP_RK3(u)::Array{Float64}
        u1 = Euler_forward(u)
        u2 = 3 / 4 * u + 1 / 4 * Euler_forward(u1)
        u_next = 1 / 3 * u + 2 / 3 * Euler_forward(u2)
    end

    while t < t_end && !quit_loop
        if t + Δt > t_end
            Δt = t_end - t
            quit_loop = true
            λ_x, λ_y = Δt / Δx, Δt / Δy
        end
        u = SSP_RK3(u)
        t += Δt
    end
    u
end


for p = 1:6
    Nₓ = 10 * 2^p
    # Grids
    x_interface = LinRange(0, 4π, Nₓ + 1)
    x = x_interface[1:end-1] .+ 2π / Nₓ
    Δx = 4π / Nₓ
    y = x

    # Initial condition
    α, β = 1 / 3, 2 / 3
    u₀(x, y) = α + β * sin((x + y) / 2)
    u0 = zeros(Float64, (Nₓ, Nₓ))
    for i = 1:Nₓ
        for j = 1:Nₓ
            # u0[i, j] = α + β * sin(x[j])
            u0[i, j] = u₀(x[j], y[i])
        end
    end

    # Solve equation
    @time u = solve(u0)
    # @profilehtml solve(u0)

    # Validation
    u_exact = zeros(Float64, (Nₓ, Nₓ))
    for i = 1:Nₓ
        for j = 1:Nₓ
            # u_exact[i, j] = exact_solver(x[j], t_end, α, β)[1]
            u_exact[i, j] = exact_solver((x[j] + y[i]) / 2, t_end, α, β)[1]
        end
    end

    path = "misc/output/t=$(t_end)/"
    mkpath(path)
    writedlm(path * "u_N=$(Nₓ).csv", u, ",")
    writedlm(path * "u_exact_N=$(Nₓ).csv", u_exact, ",")

    l1(x, y) = mean(abs.(x - y))
    l2(x, y) = sqrt(mean((x - y) .^ 2))
    l∞(x, y) = maximum(abs.(x - y))
    println(l1(u, u_exact), "\t", l2(u, u_exact), "\t", l∞(u, u_exact))
    # break
end
