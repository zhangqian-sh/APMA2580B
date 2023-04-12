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


function flux_splitting_reconstruction!(f̂⁺, f̂⁻, f⁺, f⁻)
    ϵ = 1e-6
    Nₓ = length(f⁺)
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
        f̂⁻[i+1] = WENO_u⁻([f⁺[periodic(i + k)] for k = -(stencil_size - 1):(stencil_size-1)])
        f̂⁺[i] = WENO_u⁻([f⁻[periodic(i + k)] for k = (stencil_size-1):-1:-(stencil_size - 1)])
    end
    f̂⁻[1] = f̂⁻[end]
    f̂⁺[end] = f̂⁺[1]
    f̂⁺, f̂⁻
end