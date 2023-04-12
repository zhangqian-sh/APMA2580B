using DelimitedFiles
using QuadGK
include("reconstruction.jl")
include("burgers_exact_solver.jl")

# Numerical Flux
# Here f(u) = u^2 / 2, f'(u) = u.
f(u) = u^2 / 2
# debug: check linear equation
# f(u) = u

const scheme = ARGS[1]
t_end = parse(Float64, ARGS[2])
const λ_max = 0.6 # max CFL number

# flux splitting
function LaxFriedrichsFlux(uₗ, uᵣ, α)
    (f(uₗ) + f(uᵣ) - α * (uᵣ - uₗ)) / 2
end

function f⁺(uₗ, α)
    (f(uₗ) + α * uₗ)/2
end

function f⁻(uᵣ, α)
    (f(uᵣ) - α * uᵣ)/2
end


flux = LaxFriedrichsFlux

function solve(u)
    quit_loop = false
    t = 0
    Nₓ = length(u)
    Δx = 2π / Nₓ
    Δt = λ_max * Δx^(5 / 3)
    λ = Δt / Δx

    f̂⁺, f̂⁻ = zeros(Float64, Nₓ + 1), zeros(Float64, Nₓ + 1)
    f_hat = zeros(Float64, Nₓ + 1)

    # single forward step
    function Euler_forward(u)
        v = zeros(Float64, Nₓ)
        α = maximum(abs.(u))
        # reconstruct flux split
        flux_splitting_reconstruction!(f̂⁺, f̂⁻, f⁺.(u, Ref(α)), f⁻.(u, Ref(α)))
        # monotone flux: Lax-Friedrichs
        for i = 1:Nₓ+1
            f_hat[i] = f̂⁺[i] + f̂⁻[i]
        end
        # Euler forward step
        for i = 1:Nₓ
            v[i] = u[i] - λ * (f_hat[i+1] - f_hat[i])
        end
        v
    end

    function SSP_RK3(u)
        u1 = Euler_forward(u)
        u2 = 3 / 4 * u + 1 / 4 * Euler_forward(u1)
        u_next = 1 / 3 * u + 2 / 3 * Euler_forward(u2)
    end

    while t < t_end && !quit_loop
        if t + Δt > t_end
            Δt = t_end - t
            quit_loop = true
            λ = Δt / Δx
        end
        u = SSP_RK3(u)
        t += Δt
    end
    u
end


for p = 1:6
    Nₓ = 10 * 2^p
    # Grids
    x_interface = LinRange(0, 2π, Nₓ + 1)
    x = x_interface[1:end-1] .+ π / Nₓ
    Δx = 2π / Nₓ

    # Initial condition
    α, β = 1 / 3, 2 / 3
    u₀(x) = α + β * sin(x)
    u0 = u₀.(x)
    # Project initial values to piecewise constant.
    # U₀(x) = α * x - β * cos(x) # anti-derivative
    # u0 = (U₀.(x_interface[2:end]) - U₀.(x_interface[1:end-1])) / Δx

    # Solve equation
    @time u_mean = solve(u0)

    # Validation
    # debug: check linear equation
    # exact_solver(x, t, α, β) = α + β * sin(x - t)
    # u_exact = zeros(Float64, size(u_mean))
    # for i = 1:Nₓ
    #     u_exact_mean, _err = quadgk(x -> exact_solver(x, t_end, α, β), x_interface[i], x_interface[i+1], rtol=1e-12)
    #     u_exact[i] = u_exact_mean[1] / Δx
    # end
    u_exact = exact_solver(x, t_end, α, β)

    path = "misc/output/$(scheme)/t=$(t_end)/"
    mkpath(path)
    writedlm(path * "u_N=$(Nₓ).csv", u_mean, ",")
    writedlm(path * "u_exact_N=$(Nₓ).csv", u_exact, ",")
    # break
end
