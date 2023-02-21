using DelimitedFiles
include("burgers_exact_solver.jl")

# Numerical Flux
# Here f(u) = u^2 / 2, f'(u) = u.
f(u) = u^2 / 2

function RoeFlux(uₗ, uᵣ)
    s = (uₗ + uᵣ) / 2
    if s >= 0
        f(uₗ)
    else
        f(uᵣ)
    end
end

function GodunovFlux(uₗ, uᵣ)
    if uₗ <= uᵣ
        if uₗ >= 0
            f(uₗ)
        elseif uᵣ <= 0
            f(uᵣ)
        else
            0
        end
    else
        if uᵣ >= 0
            f(uₗ)
        elseif uₗ <= 0
            f(uᵣ)
        else
            max(f(uₗ), f(uᵣ))
        end
    end
end

function LaxFriedrichsFlux(uₗ, uᵣ, u)
    α = maximum(abs.(u))
    (f(uₗ) + f(uᵣ) - α * (uᵣ - uₗ)) / 2
end

function LaxWendroffFlux(uₗ, uᵣ, λ)
    (f(uₗ) + f(uᵣ) - λ * ((uₗ + uᵣ) / 2) * (f(uᵣ) - f(uₗ))) / 2
end

α, β = 1 / 3, 2 / 3
u₀(x) = α + β * sin(x)

Nₓ = parse(Int, ARGS[2])
const Δx = 2π / Nₓ

t = parse(Float64, ARGS[1])
const λ_max = 0.6 # max CFL number
const Nₜ = ceil(t / (λ_max * 2π / Nₓ))
const Δt = t / Nₜ
const λ = Δt / Δx

flux = "LaxWendroff"

function solve(u, Nₜ)
    f = LaxWendroffFlux
    for nₜ = 1:Nₜ
        v = zeros(Float64, size(u))
        for i in eachindex(u)
            if i == 1 || i == length(u)
                fₚ, fₙ = f(u[1], u[2], λ), f(u[end-1], u[end], λ)
            else
                fₚ, fₙ = f(u[i], u[i+1], λ), f(u[i-1], u[i], λ)
            end
            v[i] = u[i] - λ * (fₚ - fₙ)
        end
        u = v
     end
    u
end

x = LinRange(0, 2π, Nₓ + 1)
u0 = u₀.(x)
u0[end] = u0[1]
@time u = solve(u0, Nₜ)

u_exact = exact_solver(x, t, α, β)

path = "misc/output/$(flux)/t=$(t)/"
mkpath(path)
writedlm(path * "u_N=$(Nₓ).csv", u, ",")
writedlm(path * "u_exact_N=$(Nₓ).csv", u_exact, ",")
