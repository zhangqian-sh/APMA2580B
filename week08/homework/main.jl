using DelimitedFiles
using QuadGK
include("burgers_exact_solver.jl")

# Numerical Flux
# Here f(u) = u^2 / 2, f'(u) = u.
f(u) = u^2 / 2

t_end = parse(Float64, ARGS[1])
const M = parse(Float64, ARGS[2])
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

function modified_minimod(a₁, a₂, a₃, C)
    if abs(a₁) <= C  # C = M * Δx^2
        a₁
    else
        minimod(a₁, a₂, a₃)
    end
end

function solve(u)
    quit_loop = false
    t = 0
    λ = λ_max
    Nₓ = length(u)
    Δx = 2π / Nₓ
    Δt = λ_max * Δx

    u⁺, u⁻ = zeros(Float64, Nₓ + 1), zeros(Float64, Nₓ + 1)

    function periodic(j)
        if 1 <= j <= Nₓ
            j
        elseif j > Nₓ
            mod(j, Nₓ)
        else
            Nₓ + j
        end
    end

    function modified_minimod_limiter(u, u⁺, u⁻)
        u_tilde_l, u_tilde_r = u - u⁺[1:end-1], u⁻[2:end] - u
        u_tilde_l_mod, u_tilde_r_mod = zeros(Float64, Nₓ), zeros(Float64, Nₓ)
        for j = 1:Nₓ
            u_tilde_l_mod[j] = modified_minimod(u_tilde_l[j], u[periodic(j + 1)] - u[j], u[j] - u[periodic(j - 1)], M * Δx^2)
            u_tilde_r_mod[j] = modified_minimod(u_tilde_r[j], u[periodic(j + 1)] - u[j], u[j] - u[periodic(j - 1)], M * Δx^2)
        end
        u⁺_mod, u⁻_mod = zeros(Float64, Nₓ + 1), zeros(Float64, Nₓ + 1)
        for j = 1:Nₓ+1
            u⁺_mod[j] = u[periodic(j)] - u_tilde_l_mod[periodic(j)]
            u⁻_mod[j] = u[periodic(j - 1)] + u_tilde_r_mod[periodic(j - 1)]
        end
        u⁺_mod, u⁻_mod
    end

    limiter = modified_minimod_limiter

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


for p = 1:6
    Nₓ = 10 * 2^p
    # Grids
    x_interface = LinRange(0, 2π, Nₓ + 1)
    x = x_interface[1:end-1] .+ π / Nₓ
    Δx = 2π / Nₓ

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

    path = "misc/output/M=$(M)/t=$(t_end)/"
    mkpath(path)
    writedlm(path * "u_N=$(Nₓ).csv", u, ",")
    writedlm(path * "u_exact_N=$(Nₓ).csv", u_exact, ",")
end
