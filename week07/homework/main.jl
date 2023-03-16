using DelimitedFiles
using QuadGK
include("burgers_exact_solver.jl")

# Numerical Flux
# Here f(u) = u^2 / 2, f'(u) = u.
f(u) = u^2 / 2

t_end = parse(Float64, ARGS[1])
Nₓ = parse(Int, ARGS[2])

const Δx = 2π / Nₓ
const λ_max = 0.6 # max CFL number

function LaxFriedrichsFlux(uₗ, uᵣ, α)
    (f(uₗ) + f(uᵣ) - α * (uᵣ - uₗ)) / 2
end

flux = LaxFriedrichsFlux

function minimod(a₁, a₂, a₃)
    if sign(a₁) == sign(a₂) == sign(a₃)
        s = sign(a₁)
        s * minimum([abs(a₁), abs(a₂), abs(a₃)])
    else
        0
    end
end

function minimod_limiter(u, u⁺, u⁻)
    u_tilde_l, u_tilde_r = u - u⁺[1:end-1], u⁻[2:end] - u
    u_tilde_l_mod, u_tilde_r_mod = zeros(Float64, Nₓ), zeros(Float64, Nₓ)
    for j = 2:Nₓ-1
        u_tilde_l_mod[j] = minimod(u_tilde_l[j], u[j+1] - u[j], u[j] - u[j-1])
        u_tilde_r_mod[j] = minimod(u_tilde_r[j], u[j+1] - u[j], u[j] - u[j-1])
    end
    u_tilde_l_mod[1] = minimod(u_tilde_l[1], u[2] - u[1], u[1] - u[end])
    u_tilde_r_mod[1] = minimod(u_tilde_r[1], u[2] - u[1], u[1] - u[end])
    u_tilde_l_mod[end] = minimod(u_tilde_l[end], u[1] - u[end], u[end] - u[end-1])
    u_tilde_r_mod[end] = minimod(u_tilde_r[end], u[1] - u[end], u[end] - u[end-1])
    u⁺_mod, u⁻_mod = zeros(Float64, Nₓ + 1), zeros(Float64, Nₓ + 1)
    for j = 2:Nₓ
        u⁺_mod[j] = u[j] - u_tilde_l_mod[j]
        u⁻_mod[j] = u[j-1] + u_tilde_r_mod[j-1]
    end
    u⁺_mod[1] = u⁺_mod[end] = u[1] - u_tilde_l_mod[1]
    u⁻_mod[1] = u⁻_mod[end] = u[end] + u_tilde_l_mod[end]
    u⁺_mod, u⁻_mod
end

limiter = minimod_limiter

function solve(u)
    quit_loop = false
    t = 0
    Δt = λ_max * Δx
    λ = λ_max
    u⁺, u⁻ = zeros(Float64, Nₓ + 1), zeros(Float64, Nₓ + 1)

    function Euler_forward(u)
        v = zeros(Float64, Nₓ)
        α = maximum(abs.(u))
        # reconstruction: stencil size = 3
        for i = 2:Nₓ-1
            u⁺[i] = 1 / 3 * u[i-1] + 5 / 6 * u[i] - 1 / 6 * u[i+1]
        end
        u⁺[Nₓ] = 1 / 3 * u[Nₓ-1] + 5 / 6 * u[Nₓ] - 1 / 6 * u[1]
        u⁺[1] = u⁺[Nₓ+1] = 1 / 3 * u[Nₓ] + 5 / 6 * u[1] - 1 / 6 * u[2] # Periodic BC
        for i = 3:Nₓ
            u⁻[i] = -1 / 6 * u[i-2] + 5 / 6 * u[i-1] + 1 / 3 * u[i]
        end
        u⁻[2] = -1 / 6 * u[Nₓ] + 5 / 6 * u[1] + 1 / 3 * u[2]
        u⁻[1] = u⁻[Nₓ+1] = -1 / 6 * u[Nₓ-1] + 5 / 6 * u[Nₓ] + 1 / 3 * u[1] # Periodic BC
        # apply limiter
        if limiter !== nothing
            u⁺, u⁻ = limiter(u, u⁺, u⁻)
        end
        # monotone flux: Lax-Friedrichs
        f_hat = zeros(Float64, Nₓ + 1)
        for i = 1:Nₓ+1
            f_hat[i] = flux(u⁻[i], u⁺[i], α)
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

# Grids
x_interface = LinRange(0, 2π, Nₓ + 1)
x = x_interface[1:end-1] .+ π / Nₓ

# Initial condition
α, β = 1 / 3, 2 / 3
u₀(x) = α + β * sin(x)
# Project initial values to piecewise constant.
U₀(x) = α * x - β * cos(x) # anti-derivative
u0 = (U₀.(x_interface[2:end]) - U₀.(x_interface[1:end-1])) / Δx

# Solve equation
@time u_mean = solve(u0)

# Reconstruct solution values on the midpoints.
u = zeros(Float64, Nₓ)
for i = 2:Nₓ-1
    u[i] = -1 / 24 * u_mean[i-1] + 26 / 24 * u_mean[i] - 1 / 24 * u_mean[i+1]
end
u[1] = -1 / 24 * u_mean[end] + 26 / 24 * u_mean[1] - 1 / 24 * u_mean[2]
u[end] = -1 / 24 * u_mean[end-1] + 26 / 24 * u_mean[end] - 1 / 24 * u_mean[1]

# Validation
u_exact = exact_solver(x, t_end, α, β)

path = "misc/output/$(limiter)/t=$(t_end)/"
mkpath(path)
writedlm(path * "u_N=$(Nₓ).csv", u, ",")
writedlm(path * "u_exact_N=$(Nₓ).csv", u_exact, ",")
