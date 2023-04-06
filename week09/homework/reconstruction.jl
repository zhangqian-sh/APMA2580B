using LinearAlgebra

function periodic_bc(j, Nₓ)
    if 1 <= j <= Nₓ
        j
    elseif j > Nₓ
        mod(j, Nₓ)
    else
        Nₓ + j
    end
end


function reconstruction_fvm!(u⁺, u⁻, u)
    Nₓ = length(u)
    periodic(j) = periodic_bc(j, Nₓ)
    # stencil size = 3
    for i = 1:Nₓ+1
        u⁺[i] = 1 / 3 * u[periodic(i - 1)] + 5 / 6 * u[periodic(i)] - 1 / 6 * u[periodic(i + 1)]
    end
    for i = 1:Nₓ+1
        u⁻[i] = -1 / 6 * u[periodic(i - 2)] + 5 / 6 * u[periodic(i - 1)] + 1 / 3 * u[periodic(i)]
    end
    u⁺, u⁻
end

function reconstruction_ENO!(u⁺, u⁻, u)
    Nₓ = length(u)
    periodic(j) = periodic_bc(j, Nₓ)
    # stencil size = 3, r = 2
    r = 2
    shift = r
    # compute undivided difference
    f = zeros(Float64, (r + 1, Nₓ + 2r))
    # initialize: 0-th order undivided difference is the value of u
    for i = -r+1:Nₓ+r
        f[1, i+shift] = u[periodic(i)]
    end
    # from low order to high order
    for k = 1:r
        for i = -r+1:Nₓ+r-k
            f[k+1, i+shift] = f[k, i+1+shift] - f[k, i+shift]
        end
    end
    # choose stencil from undivided difference
    stencil_left = Array(1:Nₓ)
    for k = 2:r+1
        for i = 1:Nₓ
            if abs(f[k, stencil_left[i]-1+shift]) < abs(f[k, stencil_left[i]+shift])
                stencil_left[i] -= 1
            end
        end
    end
    # reconstruction from the chosen stencil
    for i = 1:Nₓ
        if stencil_left[i] == i - 2
            u⁺[i] = -1 / 6 * u[periodic(i - 2)] + 5 / 6 * u[periodic(i - 1)] + 1 / 3 * u[periodic(i)]
        elseif stencil_left[i] == i - 1
            u⁺[i] = 1 / 3 * u[periodic(i - 1)] + 5 / 6 * u[periodic(i)] - 1 / 6 * u[periodic(i + 1)]
        else
            u⁺[i] = 11 / 6 * u[periodic(i)] - 7 / 6 * u[periodic(i + 1)] + 1 / 3 * u[periodic(i + 2)]
        end
    end
    u⁺[end] = u⁺[1]
    for i = 1:Nₓ
        if stencil_left[i] == i - 2
            u⁻[i+1] = 1 / 3 * u[periodic(i - 2)] - 7 / 6 * u[periodic(i - 1)] + 11 / 6 * u[periodic(i)]
        elseif stencil_left[i] == i - 1
            u⁻[i+1] = -1 / 6 * u[periodic(i - 1)] + 5 / 6 * u[periodic(i)] + 1 / 3 * u[periodic(i + 1)]
        else
            u⁻[i+1] = 1 / 3 * u[periodic(i)] + 5 / 6 * u[periodic(i + 1)] - 1 / 6 * u[periodic(i + 2)]
        end
    end
    u⁻[1] = u⁻[end]
    u⁺, u⁻
end

function reconstruction_WENO!(u⁺, u⁻, u)
    ϵ = 1e-6
    Nₓ = length(u)
    periodic(j) = periodic_bc(j, Nₓ)
    # stencil size = 3
    stencil_size = 3
    # reconstruction coefficents
    rc = [
        1/3 -7/6 11/6
        -1/6 5/6 1/3
        1/3 5/6 -1/6
    ]
    γ = [1 / 10, 3 / 5, 3 / 10]
    β_mul_symmetric = 13 / 12
    β_mul_asymmetric = 1 / 4
    β_comp_symmetric = [
        1 -2 1
        1 -2 1
        1 -2 1
    ]
    β_comp_asymmetric = [
        1 -4 3
        1 0 -1
        3 -4 1
    ]

    # temp variables
    u_temp = zeros(Float64, stencil_size)
    β = zeros(Float64, stencil_size)
    w̃ = zeros(Float64, stencil_size)

    function WENO_u⁻(u)
        for l = 1:stencil_size
            u_temp[l] = dot(rc[l, :], u[l:l+stencil_size-1])
            β[l] = β_mul_symmetric * dot(β_comp_symmetric[l, :], u[l:l+stencil_size-1])^2 +
                   β_mul_asymmetric * dot(β_comp_asymmetric[l, :], u[l:l+stencil_size-1])^2
        end
        for l = 1:stencil_size
            w̃[l] = γ[l] / (β[l] + ϵ)^2
        end
        w = w̃ / sum(w̃)
        dot(w, u_temp)
    end

    for i = 1:Nₓ
        u⁻[i+1] = WENO_u⁻([u[periodic(i + k)] for k = -(stencil_size - 1):(stencil_size-1)])
        u⁺[i] = WENO_u⁻([u[periodic(i + k)] for k = (stencil_size-1):-1:-(stencil_size - 1)])
    end
    u⁻[1] = u⁻[end]
    u⁺[end] = u⁺[1]
    u⁺, u⁻
end

# limiter
# TVD limiter
function minimod(a₁, a₂, a₃)
    if sign(a₁) == sign(a₂) == sign(a₃)
        s = sign(a₁)
        s * minimum([abs(a₁), abs(a₂), abs(a₃)])
    else
        0
    end
end

# TVB limiter
function modified_minimod(a₁, a₂, a₃, C)
    if abs(a₁) <= C  # C = M * Δx^2
        a₁
    else
        minimod(a₁, a₂, a₃)
    end
end

function modified_minimod_limiter(u, u⁺, u⁻, M)
    Nₓ = length(u)
    periodic(j) = periodic_bc(j, Nₓ)
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