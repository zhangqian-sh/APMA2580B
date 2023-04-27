include("euler_solver.jl")

function density_wave(region::ComputationalRegion)
    L = region.L
    u₀, u_exact = zeros(Float64, (5, region.N_x)), zeros(Float64, (5, region.N_x))
    ρ(x) = 1 + 0.2sin(x)

    x_interface = LinRange(0, L, region.N_x + 1)
    x = x_interface[1:end-1] .+ L / 2region.N_x

    function initial_condition()
        u₀[1, :] = ρ.(x)
        u₀[4, :] .= 1
        u₀[5, :] .= 1
        u₀
    end
    function exact_solution(t)
        u_exact[1, :] = ρ.(x .- t)
        u_exact[4, :] .= 1
        u_exact[5, :] .= 1
        u_exact
    end

    initial_condition, exact_solution
end

function burgers_like(region::ComputationalRegion)
    L = region.L
    u₀, u_exact = zeros(Float64, (5, region.N_x)), zeros(Float64, (5, region.N_x))

    Δx = L / region.N_x
    x_interface = LinRange(-Δx / 2, L - Δx / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ Δx / 2

    function initial_condition()
        u₀[1, :] = 1 .+ 0.5sin.(x)
        u₀[4, :] = 0.5sin.(x)
        u₀[5, :] = (1 .+ 0.5sin.(x)) .^ 3
        u₀
    end
    function exact_solution(t)
        w = exact_solver(x, t, √3, (1 / 2 + √3 / 2))
        z = exact_solver(x, t, -√3, (1 / 2 - √3 / 2))
        u_exact[1, :] = (w - z) / 2√3
        u_exact[4, :] = (w + z) / 2
        u_exact
    end

    initial_condition, exact_solution
end

function burgers_like_2(region::ComputationalRegion)
    L = region.L
    u₀, u_exact = zeros(Float64, (5, region.N_x)), zeros(Float64, (5, region.N_x))

    Δx = L / region.N_x
    x_interface = LinRange(-Δx / 2, L - Δx / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ Δx / 2

    function initial_condition()
        u₀[1, :] = (1 .+ 0.2sin.(x)) ./ 2√3
        u₀[4, :] = √3 .* u₀[1, :]
        u₀[5, :] = u₀[1, :] .^ 3
        u₀
    end
    function exact_solution(t)
        w = exact_solver(x, t, 1, 0.2)
        u_exact[1, :] = w ./ 2√3
        u_exact
    end

    initial_condition, exact_solution
end

function shock_tube_Sod(region::ComputationalRegion)
    L = region.L
    u₀, u_exact = zeros(Float64, (5, region.N_x)), zeros(Float64, (5, region.N_x))

    x_interface = LinRange(-L / 2, L / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ L / 2region.N_x

    function initial_condition()
        for j = 1:region.N_x
            if x[j] <= 0
                u₀[1, j] = 1.0
                u₀[4, j] = 0.0
                u₀[5, j] = 1.0
            else
                u₀[1, j] = 0.125
                u₀[4, j] = 0.0
                u₀[5, j] = 0.1
            end
        end
        u₀
    end
    function exact_solution(t)
        u_exact
    end

    initial_condition, exact_solution
end

function shock_tube_Lax(region::ComputationalRegion)
    L = region.L
    u₀, u_exact = zeros(Float64, (5, region.N_x)), zeros(Float64, (5, region.N_x))

    x_interface = LinRange(-L / 2, L / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ L / 2region.N_x

    function initial_condition()
        for j = 1:region.N_x
            if x[j] <= 0
                u₀[1, j] = 0.445
                u₀[4, j] = 0.698
                u₀[5, j] = 3.528
            else
                u₀[1, j] = 0.5
                u₀[4, j] = 0
                u₀[5, j] = 0.571
            end
        end
        u₀
    end
    function exact_solution(t)
        u_exact
    end

    initial_condition, exact_solution
end
