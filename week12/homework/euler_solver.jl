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
    L::Float64
    N_x::Int
    t_end::Float64
end

function periodic_bc(u::Array{Float64,2}, j::Int, N_x::Int)::Vector{Float64}
    if 1 <= j <= N_x
        u[:, j]
    elseif j > N_x
        u[:, mod(j, N_x)]
    else
        u[:, N_x+j]
    end
end

function shocktube_bc(u::Array{Float64,2}, j::Int, N_x::Int, uₗ::Vector{Float64}, uᵣ::Vector{Float64})::Vector{Float64}
    if 1 <= j <= N_x
        u[:, j]
    elseif j > N_x
        uᵣ[1:3]
    else
        uₗ[1:3]
    end
end

function get_boundary_condition(region::ComputationalRegion)
    if region.boundary_condition.type == "periodic"
        fn_periodic_bc(u, j) = periodic_bc(u, j, region.N_x)
    elseif region.boundary_condition.type == "shocktube"
        fn_shocktube_bc(u, j) = shocktube_bc(u, j, region.N_x, region.boundary_condition.uₗ, region.boundary_condition.uᵣ)
    else
        fn_bc(u, j) = periodic_bc(u, j, region.N_x)
    end
end

function solve(system::EulerSystem1D, region::ComputationalRegion, u::Array{T,2})::Array{T,2} where {T<:Float64}
    # load parameters
    γ = system.γ
    t_end = region.t_end
    N_x = region.N_x
    Δx = region.L / N_x
    stencil_size = 5
    half_stencil_size = Int((stencil_size + 1) / 2)

    bc = get_boundary_condition(region)

    CFL_max = 0.1
    Δt = 0

    # temp variables
    λ = ones(Float64, (3, N_x))
    L, R, Λ = zeros(Float64, (3, 3)), zeros(Float64, (3, 3)), zeros(Float64, (3, 3))
    # L, R = I, I
    fu = zeros(Float64, (3, stencil_size + 1))
    U, F = zeros(Float64, (3, stencil_size + 1)), zeros(Float64, (3, stencil_size + 1))
    F⁺, F⁻ = zeros(Float64, (3, stencil_size + 1)), zeros(Float64, (3, stencil_size + 1))
    F̂, F̂⁺, F̂⁻ = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
    f̂ = zeros(Float64, (3, N_x + 1))

    reconstruct = WENO5_reconstruction

    function f(u::Vector{Float64})::Vector{Float64}
        ρ, m, E = u
        v = m / ρ
        p = (γ - 1) * (E - 0.5 * ρ * v^2)
        [m, ρ * v^2 + p, v * (E + p)]
    end

    function compute_eigenvalue!(u)
        v = u[2, :] ./ u[1, :]
        p = (γ - 1) * (u[3, :] - 0.5 * u[1, :] .* v .^ 2)
        c = sqrt.(γ * p ./ u[1, :])
        λ[1, :], λ[2, :], λ[3, :] = v .- c, v, v .+ c
    end

    # function compute_eigenvalue!(u) end

    function compute_eigenvector!(u₀::Vector{Float64}, u₁::Vector{Float64})
        # Roe average
        ρ₀, m₀, E₀ = u₀
        ρ₁, m₁, E₁ = u₁
        v₀, v₁ = m₀ / ρ₀, m₁ / ρ₁
        p₀, p₁ = (γ - 1) * (E₀ - 0.5 * ρ₀ * v₀^2), (γ - 1) * (E₁ - 0.5 * ρ₁ * v₁^2)
        t₀ = √ρ₀ / (√ρ₀ + √ρ₁)
        t₁ = 1 - t₀
        v_xm = t₀ * v₀ + t₁ * v₁
        H₀, H₁ = (p₀ + E₀) / ρ₀, (p₁ + E₁) / ρ₁
        H_xm = t₀ * H₀ + t₁ * H₁
        qₘ = 0.5 * v_xm^2
        cₘ = sqrt((γ - 1) * (H_xm - qₘ))
        t₂ = v_xm * cₘ
        R = [1 1 1; v_xm-cₘ v_xm v_xm+cₘ; H_xm-t₂ qₘ H_xm+t₂]
        r_cm = 1 / cₘ
        b₁ = (γ - 1) * r_cm^2
        b₂ = qₘ * b₁
        t₀, t₁, t₂ = v_xm * r_cm, b₁ * v_xm, 0.5 * b₁
        L = [0.5*(b₂+t₀) -0.5*(t₁+r_cm) t₂; 1-b₂ t₁ -b₁; 0.5*(b₂-t₀) -0.5*(t₁-r_cm) t₂]
        # A = Diagonal([v_xm - cₘ, v_xm, v_xm + cₘ])
    end

    # function compute_eigenvector!(u₀::Vector{Float64}, u₁::Vector{Float64})
    #     # Raw
    #     ρ₀, m₀, E₀, v₀, p₀ = u₀
    #     ρ₁, m₁, E₁, v₁, p₁ = u₁
    #     ρ, E, v, p = (ρ₀ + ρ₁) / 2, (E₀ + E₁) / 2, (v₀ + v₁) / 2, (p₀ + p₁) / 2
    #     c = sqrt(γ * p / ρ)
    #     # H = ((E₀ + p₀) / ρ₀ + (E₁ + p₁) / ρ₁) / 2
    #     H = (E + p) / ρ
    #     R = [1 1 1; v-c v v+c; H-v*c v^2/2 H+v*c]
    #     b₁ = (γ - 1) / c^2
    #     b₂ = 0.5v^2 * b₁
    #     L = 0.5 * [b₂+v/c -(b₁ * v + 1 / c) b₁; 2(1-b₂) 2b₁*v -2b₁; b₂-v/c -(b₁ * v - 1 / c) b₁]
    # end

    # function compute_eigenvector!(u₀::Vector{Float64}, u₁::Vector{Float64}) end

    # single forward step
    function Euler_forward(u)
        compute_eigenvalue!(u)
        Λ = diagm([maximum(abs.(λ[1, :])), maximum(abs.(λ[2, :])), maximum(abs.(λ[3, :]))])
        for j = 0:N_x
            compute_eigenvector!(bc(u, j), bc(u, j + 1))
            for jₗ = 1:stencil_size+1
                u_l = bc(u, j + jₗ - half_stencil_size)
                U[:, jₗ] = L * u_l
                fu[:, jₗ] = f(u_l)
                F[:, jₗ] = L * fu[:, jₗ]
                # F⁺[:, jₗ] = 0.5 * (F[:, jₗ] + Λ * U[:, jₗ])
                # F⁻[:, jₗ] = 0.5 * (F[:, jₗ] - Λ * U[:, jₗ])
            end
            F⁺, F⁻ = 0.5 * (F + Λ * U), 0.5 * (F - Λ * U)
            for k = 1:3
                F̂⁻[k] = reconstruct(F⁺[k, 1], F⁺[k, 2], F⁺[k, 3], F⁺[k, 4], F⁺[k, 5])
                F̂⁺[k] = reconstruct(F⁻[k, 6], F⁻[k, 5], F⁻[k, 4], F⁻[k, 3], F⁻[k, 2])
            end
            F̂ = F̂⁺ + F̂⁻
            f̂[:, j+1] = R * F̂
        end

        # u_euler_forward = zeros(Float64, size(u))
        # u_euler_forward = u - Δt / Δx * (f̂[:, 2:end] - f̂[:, 1:end-1])
        # u_euler_forward
        -(f̂[:, 2:end] - f̂[:, 1:end-1]) / Δx
    end

    # function SSP_RK3(u)
    #     u1 = Euler_forward(u)
    #     u2 = 3 / 4 * u + 1 / 4 * Euler_forward(u1)
    #     u_next = 1 / 3 * u + 2 / 3 * Euler_forward(u2)
    # end

    # f1, f2, f3, f4 = zeros(Float64, (3, N_x)), zeros(Float64, (3, N_x)), zeros(Float64, (3, N_x)), zeros(Float64, (3, N_x))

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

    function SSP_RK4(u)
        u1 = u + 0.39175222700392 * Δt * Euler_forward(u)
        u2 = 0.44437049406734 * u + 0.55562950593266 * u1 +
             0.36841059262959 * Δt * Euler_forward(u1)
        u3 = 0.62010185138540 * u + 0.37989814861460 * u2 +
             0.25189177424738 * Δt * Euler_forward(u2)
        u4 = 0.17807995410773 * u + 0.82192004589227 * u3 +
             0.54497475021237 * Δt * Euler_forward(u3)
        w = 0.00683325884039 * u + 0.51723167208978 * u2 +
            0.12759831133288 * u3 + 0.34833675773694 * u4 +
            0.08460416338212 * Δt * Euler_forward(u3) +
            0.22600748319395 * Δt * Euler_forward(u4)
    end

    t = 0
    quit_loop = false

    while t < t_end && !quit_loop
        # determine CFL number
        compute_eigenvalue!(u)
        α = maximum(abs.(λ))
        CFL = CFL_max / α
        Δt = CFL_max * (Δx)^(1)
        # Δt = 0.1*Δx
        # println(t)
        if t + Δt > t_end
            Δt = t_end - t
            quit_loop = true
        end
        # u = RK4(u)
        u = SSP_RK3(u)
        t += Δt
    end
    u
end