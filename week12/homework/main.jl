using DelimitedFiles
using QuadGK
using Statistics
using StatProfilerHTML

include("data.jl")
include("euler_solver.jl")
include("burgers_exact_solver.jl")

data = ARGS[1]
t_end = parse(Float64, ARGS[2])

function compute_non_conservative_variable!(γ, u)
    u[4, :] = u[2, :] ./ u[1, :]
    u[5, :] = (γ - 1) * (u[3, :] .- 0.5 * u[1, :] .* u[4, :] .^ 2)
end

function compute_conservative_variable!(γ, u)
    u[2, :] = u[1, :] .* u[4, :]
    u[3, :] = u[5, :] / (γ - 1) + 0.5 * u[1, :] .* u[4, :] .^ 2
end


for p = 1:6
    N_x = 10 * 2^p

    if data == "density_wave"
        system = EulerSystem1D(1.4)
        region = ComputationalRegion(PeriodicBoundaryCondition("periodic"), 2π, N_x, t_end)
        initial_condition, exact_solution = density_wave(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "burgers_like"
        system = EulerSystem1D(3)
        region = ComputationalRegion(PeriodicBoundaryCondition("periodic"), 2π, N_x, t_end)
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
        uₗ, uᵣ = [1.0, 0, 0, 0, 1.0], [0.125, 0, 0, 0, 0.1]
        compute_conservative_variable!(system.γ, uₗ)
        compute_conservative_variable!(system.γ, uᵣ)
        region = ComputationalRegion(ShockTubeBoundaryCondition("shocktube", uₗ, uᵣ), 10.0, N_x, t_end)
        initial_condition, exact_solution = shock_tube_Sod(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    elseif data == "Lax"
        system = EulerSystem1D(1.4)
        uₗ, uᵣ = [0.445, 0, 0, 0.698, 3.528], [0.5, 0, 0, 0, 0.571]
        compute_conservative_variable!(system.γ, uₗ)
        compute_conservative_variable!(system.γ, uᵣ)
        region = ComputationalRegion(ShockTubeBoundaryCondition("shocktube", uₗ, uᵣ), 10.0, N_x, t_end)
        initial_condition, exact_solution = shock_tube_Lax(region)
        u0 = initial_condition()
        compute_conservative_variable!(system.γ, u0)
    end

    # Solve equation
    u = zeros(Float64, size(u0))
    u0 = u0[1:3, :]
    @time u[1:3, :] = solve(system, region, u0)
    compute_non_conservative_variable!(system.γ, u)
    # @profilehtml solve(system, region, u0)

    u_exact = exact_solution(t_end)
    compute_conservative_variable!(system.γ, u_exact)

    path = "misc/output/$(data)/t=$(t_end)/"
    mkpath(path)
    writedlm(path * "u_N=$(N_x).csv", u, ",")
    writedlm(path * "u_exact_N=$(N_x).csv", u_exact, ",")
    l1(x, y) = mean(abs.(x - y))
    l2(x, y) = sqrt(mean((x - y) .^ 2))
    l∞(x, y) = maximum(abs.(x - y))
    println(l1(u, u_exact), "\t", l2(u, u_exact), "\t", l∞(u, u_exact))
    if data == "burgers_like"
        println(l1(u[[1, 4], :], u_exact[[1, 4], :]), "\t", l2(u[[1, 4], :], u_exact[[1, 4], :]), "\t", l∞(u[[1, 4], :], u_exact[[1, 4], :]))
    elseif data == "burgers_like_2"
        println(l1(u[1, :], u_exact[1, :]), "\t", l2(u[1, :], u_exact[1, :]), "\t", l∞(u[1, :], u_exact[1, :]))
    end
    # break
    flush(stdout)
end
