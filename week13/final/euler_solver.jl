include("reconstruction.jl")

struct EulerSystem1D
    γ::Float64
end

abstract type BoundaryCondition end

struct PeriodicBoundaryCondition <: BoundaryCondition
    type::String
end

struct ShockTubeBoundaryCondition <: BoundaryCondition
    type::String
    uₗ::Vector{Float64}
    uᵣ::Vector{Float64}
end

struct ComputationalRegion
    boundary_condition::BoundaryCondition
    L_x::Float64
    L_y::Float64
    N_x::Int
    N_y::Int
    t_end::Float64
end

function periodic_bc(u::Array{Float64,3}, i::Int, j::Int, N_y::Int, N_x::Int)::Vector{Float64}
    if 1 <= j <= N_x
        idx_j = j
    elseif j > N_x
        idx_j = mod(j, N_x)
    else
        idx_j = N_x + j
    end
    if 1 <= i <= N_y
        idx_i = i
    elseif i > N_y
        idx_i = mod(i, N_y)
    else
        idx_i = N_y + i
    end
    u[:, idx_i, idx_j]
end

function shocktube_bc(u::Array{Float64,3}, i::Int, j::Int, N_y::Int, N_x::Int, uₗ::Vector{Float64}, uᵣ::Vector{Float64})::Vector{Float64}
    # println(i, " ", j)
    idx_i, idx_j = 0, 0
    if i + j > N_x + N_y
        return uᵣ[1:4]
    elseif i + j < 2
        return uₗ[1:4]
    else
        idx_i, idx_j = i, j
    end
    # println(idx_i, " ", idx_j)
    if 1 <= idx_j <= N_x && 1 <= idx_i <= N_y
        idx_i, idx_j = idx_i, idx_j
    elseif idx_j > N_x || idx_i < 1 || idx_i > N_y || idx_j < 1
        idx_i, idx_j = Int(ceil((idx_i + idx_j) / 2)), Int(floor((idx_i + idx_j) / 2))
    else
        throw(ErrorException("i, j are not in the correct domain"))
    end
    u[:, idx_i, idx_j]
end

function get_boundary_condition(region::ComputationalRegion)
    if region.boundary_condition.type == "periodic"
        fn_periodic_bc(u, i, j) = periodic_bc(u, i, j, region.N_y, region.N_x)
    elseif region.boundary_condition.type == "shocktube"
        fn_shocktube_bc(u, i, j) = shocktube_bc(u, i, j, region.N_y, region.N_x, region.boundary_condition.uₗ, region.boundary_condition.uᵣ)
    else
        fn_bc(u, j) = periodic_bc(u, i, j, region.N_y, region.N_x)
    end
end

function solve(system::EulerSystem1D, region::ComputationalRegion, u::Array{Float64,3})::Array{Float64,3}
    # load parameters
    γ = system.γ
    t_end = region.t_end
    N_x, N_y = region.N_x, region.N_y
    Δx, Δy = region.L_x / N_x, region.L_y / N_y
    stencil_size = 5
    half_stencil_size = Int((stencil_size + 1) / 2)

    bc = get_boundary_condition(region)

    CFL_max = 0.4
    Δt = 0

    # temp variables
    λ_f, λ_g = ones(Float64, (4, N_x)), ones(Float64, (4, N_y))
    λ = ones(Float64, (2, 4, N_y, N_x))
    # L, R, Λ = zeros(Float64, (4, 4)), zeros(Float64, (4, 4)), zeros(Float64, (4, 4))
    L, R = I, I
    U = zeros(Float64, (4, stencil_size + 1))
    F, F⁺, F⁻ = zeros(Float64, (4, stencil_size + 1)), zeros(Float64, (4, stencil_size + 1)), zeros(Float64, (4, stencil_size + 1))
    F̂, F̂⁺, F̂⁻ = zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4)
    G, G⁺, G⁻ = zeros(Float64, (4, stencil_size + 1)), zeros(Float64, (4, stencil_size + 1)), zeros(Float64, (4, stencil_size + 1))
    Ĝ, Ĝ⁺, Ĝ⁻ = zeros(Float64, 4), zeros(Float64, 4), zeros(Float64, 4)
    f̂, ĝ = zeros(Float64, (4, N_y, N_x + 1)), zeros(Float64, (4, N_y + 1, N_x))

    reconstruct = WENO5_reconstruction

    function f(u::Vector{Float64})::Vector{Float64}
        ρ, ρv, ρw, E = u
        v, w = ρv / ρ, ρw / ρ
        p = (γ - 1) * (E - 0.5 * ρ * (v^2 + w^2))
        [ρv, ρ * v^2 + p, ρ * v * w, v * (E + p)]
    end

    function g(u::Vector{Float64})::Vector{Float64}
        ρ, ρv, ρw, E = u
        v, w = ρv / ρ, ρw / ρ
        p = (γ - 1) * (E - 0.5 * ρ * (v^2 + w^2))
        [ρw, ρ * v * w, ρ * w^2 + p, w * (E + p)]
    end

    function compute_eigenvalue!(u::Array{Float64,3})
        ρ, ρv, ρw, E = u[1, :, :], u[2, :, :], u[3, :, :], u[4, :, :]
        v, w = ρv ./ ρ, ρw ./ ρ
        p = @. (γ - 1) * (E - 0.5 * ρ * (v^2 + w^2))
        c = @. sqrt(γ * p / ρ)
        λ[1, 1, :, :], λ[1, 2, :, :], λ[1, 3, :, :], λ[1, 4, :, :] = v - c, v, v, v + c
        λ[2, 1, :, :], λ[2, 2, :, :], λ[2, 3, :, :], λ[2, 4, :, :] = w - c, w, w, w + c
    end

    function compute_f_eigenvalue!(u::Array{Float64,2})
        ρ, ρv, ρw, E = u[1, :], u[2, :], u[3, :], u[4, :]
        v, w = ρv ./ ρ, ρw ./ ρ
        p = @. (γ - 1) * (E - 0.5 * ρ * (v^2 + w^2))
        c = @. sqrt(γ * p / ρ)
        λ_f[1, :], λ_f[2, :], λ_f[3, :], λ_f[4, :] = v - c, v, v, v + c
    end

    function compute_f_eigenvector!(u₀::Vector{Float64}, u₁::Vector{Float64})
        ρ₀, ρv₀, ρw₀, E₀ = u₀
        ρ₁, ρv₁, ρw₁, E₁ = u₁
        # v₀, v₁ = ρv₀ / ρ₀, ρv₁ / ρ₁
        # w₀, w₁ = ρw₀ / ρ₀, ρw₁ / ρ₁
        # p₀, p₁ = (γ - 1) * (E₀ - 0.5 * ρ₀ * (v₀^2 + w₀^2)), (γ - 1) * (E₁ - 0.5 * ρ₁ * (v₁^2 + w₁^2))
        # c₀, c₁ = sqrt(γ * p₀ / ρ₀), sqrt(γ * p₁ / ρ₁)
        # H₀, H₁ = (p₀ + E₀) / ρ₀, (p₁ + E₁) / ρ₁
        # if debug
        #     @show ρ₀, ρv₀, ρw₀, E₀, ρ₁, ρv₁, ρw₁, E₁
        # end
        # Naive average
        ρ = (ρ₀ + ρ₁) / 2
        ρv, ρw = (ρv₀ + ρv₁) / 2, (ρw₀ + ρw₁) / 2
        v_xm, w_xm = ρv / ρ, ρw / ρ
        E = (E₀ + E₁) / 2
        p = (γ - 1) * (E - 0.5 * (ρv^2 + ρw^2) / ρ)
        cₘ = sqrt(γ * p / ρ)
        H_xm = (E + p) / ρ
        # # Roe average
        # t₀ = √ρ₀ / (√ρ₀ + √ρ₁)
        # t₁ = 1 - t₀
        # v_xm = t₀ * v₀ + t₁ * v₁
        # w_xm = t₀ * w₀ + t₁ * w₁
        # H_xm = t₀ * H₀ + t₁ * H₁
        # cₘ = sqrt((γ - 1) * (H_xm - qₘ))
        qₘ = 0.5 * (v_xm^2 + w_xm^2)
        R = [1 1 0 1
            v_xm-cₘ v_xm 0 v_xm+cₘ
            w_xm w_xm 1 w_xm
            H_xm-v_xm*cₘ qₘ w_xm H_xm+v_xm*cₘ]
        L = [(H_xm+cₘ/(γ-1)*(v_xm-cₘ))*(γ-1)/(2*cₘ^2) -(v_xm + cₘ / (γ - 1))*(γ-1)/(2cₘ^2) -w_xm*(γ-1)/(2cₘ^2) 1*(γ-1)/(2cₘ^2)
            (-2H_xm+4/(γ-1)*cₘ^2)*(γ-1)/(2cₘ^2) 2v_xm*(γ-1)/(2cₘ^2) 2w_xm*(γ-1)/(2cₘ^2) -2*(γ-1)/(2cₘ^2)
            -2w_xm*cₘ^2/(γ-1)*(γ-1)/(2cₘ^2) 0 1 0
            (H_xm-cₘ/(γ-1)*(v_xm+cₘ))*(γ-1)/(2cₘ^2) (-v_xm+cₘ/(γ-1))*(γ-1)/(2cₘ^2) -w_xm*(γ-1)/(2cₘ^2) 1*(γ-1)/(2cₘ^2)]
        # L = inv(R)
        # if debug
        #     display(L * R)
        #     display(R * L)
        # end
    end

    function compute_g_eigenvalue!(u::Array{Float64,2})
        ρ, ρv, ρw, E = u[1, :], u[2, :], u[3, :], u[4, :]
        v, w = ρv ./ ρ, ρw ./ ρ
        p = @. (γ - 1) * (E - 0.5 * ρ * (v^2 + w^2))
        c = @. sqrt(γ * p / ρ)
        λ_g[1, :], λ_g[2, :], λ_g[3, :], λ_g[4, :] = w - c, w, w, w + c
    end

    function compute_g_eigenvector!(u₀::Vector{Float64}, u₁::Vector{Float64})
        ρ₀, ρv₀, ρw₀, E₀ = u₀
        ρ₁, ρv₁, ρw₁, E₁ = u₁
        # v₀, v₁ = ρv₀ / ρ₀, ρv₁ / ρ₁
        # w₀, w₁ = ρw₀ / ρ₀, ρw₁ / ρ₁
        # p₀, p₁ = (γ - 1) * (E₀ - 0.5 * ρ₀ * (v₀^2 + w₀^2)), (γ - 1) * (E₁ - 0.5 * ρ₁ * (v₁^2 + w₁^2))
        # c₀, c₁ = sqrt(γ * p₀ / ρ₀), sqrt(γ * p₁ / ρ₁)
        # H₀, H₁ = (p₀ + E₀) / ρ₀, (p₁ + E₁) / ρ₁
        # Naive average
        ρ = (ρ₀ + ρ₁) / 2
        ρv, ρw = (ρv₀ + ρv₁) / 2, (ρw₀ + ρw₁) / 2
        v_xm, w_xm = ρv / ρ, ρw / ρ
        E = (E₀ + E₁) / 2
        p = (γ - 1) * (E - 0.5 * (ρv^2 + ρw^2) / ρ)
        cₘ = sqrt(γ * p / ρ)
        H_xm = (E + p) / ρ
        # # Roe average
        # t₀ = √ρ₀ / (√ρ₀ + √ρ₁)
        # t₁ = 1 - t₀
        # v_xm = t₀ * v₀ + t₁ * v₁
        # w_xm = t₀ * w₀ + t₁ * w₁
        # H_xm = t₀ * H₀ + t₁ * H₁
        # cₘ = sqrt((γ - 1) * (H_xm - qₘ))
        qₘ = 0.5 * (v_xm^2 + w_xm^2)
        R = [1 1 0 1
            v_xm v_xm 1 v_xm
            w_xm-cₘ w_xm 0 w_xm+cₘ
            H_xm-w_xm*cₘ qₘ v_xm H_xm+w_xm*cₘ]
        L = [(H_xm+cₘ/(γ-1)*(w_xm-cₘ))*(γ-1)/(2cₘ^2) -v_xm*(γ-1)/(2cₘ^2) -(w_xm + cₘ / (γ - 1))*(γ-1)/(2cₘ^2) 1*(γ-1)/(2cₘ^2)
            (-2H_xm+4/(γ-1)*cₘ^2)*(γ-1)/(2cₘ^2) 2v_xm*(γ-1)/(2cₘ^2) 2w_xm*(γ-1)/(2cₘ^2) -2*(γ-1)/(2cₘ^2)
            -2v_xm*cₘ^2/(γ-1)*(γ-1)/(2cₘ^2) 1 0 0
            (H_xm-cₘ/(γ-1)*(w_xm+cₘ))*(γ-1)/(2cₘ^2) -v_xm*(γ-1)/(2cₘ^2) (-w_xm+cₘ/(γ-1))*(γ-1)/(2cₘ^2) 1*(γ-1)/(2cₘ^2)]
        # L = inv(R)
    end

    # single forward step
    function Euler_forward(u)
        # println("t = ", t)
        # display(u[4, 8:12, 8:12])
        for i = 1:N_y
            compute_f_eigenvalue!(u[:, i, :])
            Λ = diagm([maximum(abs.(λ_f[1, :])), maximum(abs.(λ_f[2, :])), maximum(abs.(λ_f[3, :])), maximum(abs.(λ_f[4, :]))])
            for j = 0:N_x
                compute_f_eigenvector!(bc(u, i, j), bc(u, i, j + 1))
                for jₗ = 1:stencil_size+1
                    u_l = bc(u, i, j + jₗ - half_stencil_size)
                    U[:, jₗ] = L * u_l
                    F[:, jₗ] = L * f(u_l)
                end
                F⁺, F⁻ = 0.5 * (F + Λ * U), 0.5 * (F - Λ * U)
                for k = 1:4
                    F̂⁻[k] = reconstruct(F⁺[k, 1], F⁺[k, 2], F⁺[k, 3], F⁺[k, 4], F⁺[k, 5])
                    F̂⁺[k] = reconstruct(F⁻[k, 6], F⁻[k, 5], F⁻[k, 4], F⁻[k, 3], F⁻[k, 2])
                end
                F̂ = F̂⁺ + F̂⁻
                f̂[:, i, j+1] = R * F̂
            end
        end
        for j = 1:N_x
            compute_g_eigenvalue!(u[:, :, j])
            Λ = diagm([maximum(abs.(λ_g[1, :])), maximum(abs.(λ_g[2, :])), maximum(abs.(λ_g[3, :])), maximum(abs.(λ_g[4, :]))])
            for i = 0:N_y
                compute_g_eigenvector!(bc(u, i, j), bc(u, i + 1, j))
                for iₗ = 1:stencil_size+1
                    u_l = bc(u, i + iₗ - half_stencil_size, j)
                    U[:, iₗ] = L * u_l
                    G[:, iₗ] = L * g(u_l)
                end
                G⁺, G⁻ = 0.5 * (G + Λ * U), 0.5 * (G - Λ * U)
                for k = 1:4
                    Ĝ⁻[k] = reconstruct(G⁺[k, 1], G⁺[k, 2], G⁺[k, 3], G⁺[k, 4], G⁺[k, 5])
                    Ĝ⁺[k] = reconstruct(G⁻[k, 6], G⁻[k, 5], G⁻[k, 4], G⁻[k, 3], G⁻[k, 2])
                end
                Ĝ = Ĝ⁺ + Ĝ⁻
                ĝ[:, i+1, j] = R * Ĝ
            end
        end
        -(f̂[:, :, 2:end] - f̂[:, :, 1:end-1]) / Δx - (ĝ[:, 2:end, :] - ĝ[:, 1:end-1, :]) / Δy
    end

    function SSP_RK3(u)
        f1 = Δt * Euler_forward(u)
        f2 = Δt * Euler_forward(u + f1)
        f3 = Δt * Euler_forward(u + 1 / 4 * f1 + 1 / 4 * f2)
        u_next = u + 1 / 6 * f1 + 1 / 6 * f2 + 2 / 3 * f3
    end

    function RK4(u)
        f1 = Δt * Euler_forward(u)
        f2 = Δt * Euler_forward(u + 1 / 2 * f1)
        f3 = Δt * Euler_forward(u + 1 / 2 * f2)
        f4 = Δt * Euler_forward(u + f3)
        u_next = u + 1 / 6 * f1 + 1 / 3 * f2 + 1 / 3 * f3 + 1 / 6 * f4
    end

    t = 0
    quit_loop = false
    step = 0
    while t < t_end && !quit_loop
        # determine CFL number
        compute_eigenvalue!(u)
        α = maximum(abs.(λ[1, :, :, :])) + maximum(abs.(λ[2, :, :, :]))
        # @show t, α
        CFL = CFL_max / α
        Δt = CFL_max * (Δx)^(1.25)
        if step % 100 == 0
            @show t, α, Δt
            flush(stdout)
        end
        step += 1
        if t + Δt > t_end
            Δt = t_end - t
            quit_loop = true
        end
        u = RK4(u)
        # u = SSP_RK3(u)
        t += Δt
    end
    u
end