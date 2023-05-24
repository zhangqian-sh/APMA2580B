include("euler_solver.jl")

function density_wave(region::ComputationalRegion)
    L_x, L_y = region.L_x, region.L_y
    u₀, u_exact = zeros(Float64, (7, region.N_y, region.N_x)), zeros(Float64, (7, region.N_y, region.N_x))
    ρ(x, y) = 1 + 0.2sin(x + y)

    # Grids
    x_interface = LinRange(0, L_x, region.N_x + 1)
    x = x_interface[1:end-1] .+ L_x / 2region.N_x
    Δx = L_x / region.N_x
    y_interface = LinRange(0, L_y, region.N_y + 1)
    y = y_interface[1:end-1] .+ L_y / 2region.N_y
    Δy = L_y / region.N_y

    function initial_condition()
        for i = 1:region.N_y
            for j = 1:region.N_x
                u₀[1, i, j] = ρ(x[j], y[i])
                u₀[5, i, j] = 1.0
                u₀[6, i, j] = 1.0
                u₀[7, i, j] = 1.0
            end
        end
        u₀
    end
    function exact_solution(t)
        for i = 1:region.N_y
            for j = 1:region.N_x
                u_exact[1, i, j] = ρ(x[j] - t, y[i] - t)
                u_exact[5, i, j] = 1.0
                u_exact[6, i, j] = 1.0
                u_exact[7, i, j] = 1.0
            end
        end
        u_exact
    end

    initial_condition, exact_solution
end

function burgers_like(region::ComputationalRegion)
    L_x, L_y = region.L_x, region.L_y
    u₀, u_exact = zeros(Float64, (7, region.N_y, region.N_x)), zeros(Float64, (7, region.N_y, region.N_x))

    # Grids
    x_interface = LinRange(0, L_x, region.N_x + 1)
    x = x_interface[1:end-1] .+ L_x / 2region.N_x
    Δx = L_x / region.N_x
    y_interface = LinRange(0, L_y, region.N_y + 1)
    y = y_interface[1:end-1] .+ L_y / 2region.N_y
    Δy = L_y / region.N_y

    ρ₀(x) = 1 + 1 / 2 * sin(x)
    v₀(x) = 1 / 2 + sin(x)
    p₀(x) = (1 + 1 / 2 * sin(x))^3

    coordinate_transform(x, y) = x / √2 + y / √2


    function initial_condition()
        for i = 1:region.N_y
            for j = 1:region.N_x
                u₀[1, i, j] = ρ₀(coordinate_transform(x[j], y[i]))
                u₀[5, i, j] = v₀(coordinate_transform(x[j], y[i])) / √2
                u₀[6, i, j] = v₀(coordinate_transform(x[j], y[i])) / √2
                u₀[7, i, j] = p₀(coordinate_transform(x[j], y[i]))
            end
        end
        u₀
    end
    function exact_solution(t)
        for i = 1:region.N_y
            for j = 1:region.N_x
                w = exact_solver(coordinate_transform(x[j], y[i]), t, 1 / 2 + √3, (1 + √3 / 2))[1]
                z = exact_solver(coordinate_transform(x[j], y[i]), t, 1 / 2 - √3, (1 - √3 / 2))[1]
                u_exact[1, i, j] = (w - z) / 2√3
                u_exact[5, i, j] = (w + z) / 2√2
                u_exact[6, i, j] = (w + z) / 2√2
            end
        end
        u_exact

        # w = exact_solver(x, t, √3, (1 / 2 + √3 / 2))
        # z = exact_solver(x, t, -√3, (1 / 2 - √3 / 2))
        # u_exact[1, :] = (w - z) / 2√3
        # u_exact[5, :] = (w + z) / 2
        # u_exact
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
    L_x, L_y = region.L_x, region.L_y
    u₀, u_exact = zeros(Float64, (7, region.N_y, region.N_x)), zeros(Float64, (7, region.N_y, region.N_x))

    # Grids
    x_interface = LinRange(-L_x / 2, L_x / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ L_x / 2region.N_x
    Δx = L_x / region.N_x
    y_interface = LinRange(-L_y / 2, L_y / 2, region.N_y + 1)
    y = y_interface[1:end-1] .+ L_y / 2region.N_y
    Δy = L_y / region.N_y

    coordinate_transform(x, y) = x / √2 + y / √2

    function initial_condition()
        for i = 1:region.N_y
            for j = 1:region.N_x
                if coordinate_transform(x[j], y[i]) <= 0
                    u₀[1, i, j] = 1.0
                    u₀[5, i, j] = 0.0
                    u₀[6, i, j] = 0.0
                    u₀[7, i, j] = 1.0
                else
                    u₀[1, i, j] = 0.125
                    u₀[5, i, j] = 0.0
                    u₀[6, i, j] = 0.0
                    u₀[7, i, j] = 0.1
                end
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
    L_x, L_y = region.L_x, region.L_y
    u₀, u_exact = zeros(Float64, (7, region.N_y, region.N_x)), zeros(Float64, (7, region.N_y, region.N_x))

    # Grids
    x_interface = LinRange(-L_x / 2, L_x / 2, region.N_x + 1)
    x = x_interface[1:end-1] .+ L_x / 2region.N_x
    Δx = L_x / region.N_x
    y_interface = LinRange(-L_y / 2, L_y / 2, region.N_y + 1)
    y = y_interface[1:end-1] .+ L_y / 2region.N_y
    Δy = L_y / region.N_y

    coordinate_transform(x, y) = x / √2 + y / √2

    function initial_condition()
        for i = 1:region.N_y
            for j = 1:region.N_x
                if coordinate_transform(x[j], y[i]) <= 0
                    u₀[1, i, j] = 0.445
                    u₀[5, i, j] = 0.698 / √2
                    u₀[6, i, j] = 0.698 / √2
                    u₀[7, i, j] = 3.528
                else
                    u₀[1, i, j] = 0.5
                    u₀[5, i, j] = 0.0
                    u₀[6, i, j] = 0.0
                    u₀[7, i, j] = 0.571
                end
            end
        end
        u₀
    end
    function exact_solution(t)
        u_exact
    end

    initial_condition, exact_solution
end
