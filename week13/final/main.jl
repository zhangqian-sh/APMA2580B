using DelimitedFiles
using QuadGK
using Statistics
using JLD2
# using StatProfilerHTML

include("data.jl")
include("euler_solver.jl")
include("burgers_exact_solver.jl")

data = ARGS[1]
t_end = parse(Float64, ARGS[2])

# u = [ρ, ρv, ρw, E, v, w, p]
function compute_non_conservative_variable!(γ, u)
    ρ, ρv, ρw, E = u[1, :, :], u[2, :, :], u[3, :, :], u[4, :, :]
    v, w = ρv ./ ρ, ρw ./ ρ
    u[5, :, :], u[6, :, :] = v, w
    u[7, :, :] = (γ - 1) * (E .- 0.5 * ρ .* (v .^ 2 + w .^ 2))
end

function compute_conservative_variable!(γ, u)
    ρ, v, w, p = u[1, :, :], u[5, :, :], u[6, :, :], u[7, :, :]
    u[2, :, :], u[3, :, :] = ρ .* v, ρ .* w
    u[4, :, :] = p / (γ - 1) + 0.5 * ρ .* (v .^ 2 + w .^ 2)
end


for p = 1:6
    println("Grid num = $p")
    N_x = 10 * 2^p
    N_y = N_x

    if data == "density_wave"
        system = EulerSystem1D(1.4)
        region = ComputationalRegion(PeriodicBoundaryCondition("periodic"), 2π, 2π, N_y, N_x, t_end)
        initial_condition, exact_solution = density_wave(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "burgers_like"
        system = EulerSystem1D(3)
        region = ComputationalRegion(PeriodicBoundaryCondition("periodic"), 2√2π, 2√2π, N_y, N_x, t_end)
        initial_condition, exact_solution = burgers_like(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "burgers_like_2"
        system = EulerSystem1D(3)
        region = ComputationalRegion(PeriodicBoundaryCondition("periodic"), 2π, N_x, t_end)
        initial_condition, exact_solution = burgers_like_2(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "Sod"
        system = EulerSystem1D(1.4)
        uₗ, uᵣ = [1.0, 0, 0, 0, 0, 0, 1.0], [0.125, 0, 0, 0, 0, 0, 0.1]
        compute_conservative_variable!(system.γ, uₗ)
        compute_conservative_variable!(system.γ, uᵣ)
        region = ComputationalRegion(ShockTubeBoundaryCondition("shocktube", uₗ, uᵣ), 5√2, 5√2, N_y, N_x, t_end)
        initial_condition, exact_solution = shock_tube_Sod(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "Lax"
        system = EulerSystem1D(1.4)
        uₗ, uᵣ = [0.445, 0, 0, 0, 0.698 / √2, 0.698 / √2, 3.528], [0.5, 0, 0, 0, 0, 0, 0.571]
        compute_conservative_variable!(system.γ, uₗ)
        compute_conservative_variable!(system.γ, uᵣ)
        region = ComputationalRegion(ShockTubeBoundaryCondition("shocktube", uₗ, uᵣ), 5√2, 5√2, N_y, N_x, t_end)
        initial_condition, exact_solution = shock_tube_Lax(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    end

    # Solve equation
    u = zeros(Float64, size(u0))
    u0 = u0[1:4, :, :]
    @time u[1:4, :, :] = solve(system, region, u0)
    compute_non_conservative_variable!(system.γ, u)
    # @profilehtml solve(system, region, u0)

    u_exact = exact_solution(t_end)
    compute_conservative_variable!(system.γ, u_exact)

    path = "misc/output/$(data)-test/t=$(t_end)/"
    mkpath(path)
    jldsave(path * "u_N=$(N_x).jld2";
        rho=u[1, :, :], rhov=u[2, :, :], rhow=u[3, :, :], E=u[4, :, :],
        v=u[5, :, :], w=u[6, :, :], p=u[7, :, :]
    )
    jldsave(path * "u_exact_N=$(N_x).jld2";
        rho=u_exact[1, :, :], rhov=u_exact[2, :, :], rhow=u_exact[3, :, :], E=u_exact[4, :, :],
        v=u_exact[5, :, :], w=u_exact[6, :, :], p=u_exact[7, :, :]
    )
    # writedlm(path * "u_N=$(N_x).csv", u[1, :, :], ",")
    # writedlm(path * "u_exact_N=$(N_x).csv", u_exact[1, :, :], ",")
    l1(x, y) = mean(abs.(x - y))
    l2(x, y) = sqrt(mean((x - y) .^ 2))
    l∞(x, y) = maximum(abs.(x - y))
    println(l1(u, u_exact), "\t", l2(u, u_exact), "\t", l∞(u, u_exact))
    if data == "burgers_like"
        println(l1(u[[1, 5, 6], :, :], u_exact[[1, 5, 6], :, :]), "\t", l2(u[[1, 5, 6], :, :], u_exact[[1, 5, 6], :, :]), "\t", l∞(u[[1, 5, 6], :, :], u_exact[[1, 5, 6], :, :]))
    elseif data == "burgers_like_2"
        println(l1(u[1, :], u_exact[1, :]), "\t", l2(u[1, :], u_exact[1, :]), "\t", l∞(u[1, :], u_exact[1, :]))
    end
    # break
    flush(stdout)
end
