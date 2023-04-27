using LinearAlgebra

function periodic_bc(j, Nₓ)::Int
    if 1 <= j <= Nₓ
        j
    elseif j > Nₓ
        mod(j, Nₓ)
    else
        Nₓ + j
    end
end

function WENO5_reconstruct(u₁, u₂, u₃, u₄, u₅)::Float64
    ϵ = 1e-6
    v₁ = 1 / 3 * u₁ - 7 / 6 * u₂ + 11 / 6 * u₃
    v₂ = -1 / 6 * u₂ + 5 / 6 * u₃ + 1 / 3 * u₄
    v₃ = 1 / 3 * u₃ + 5 / 6 * u₄ - 1 / 6 * u₅
    β₁ = 13 / 12 * (u₁ - 2u₂ + u₃)^2 + 1 / 4 * (u₁ - 4u₂ + 3u₃)^2
    β₂ = 13 / 12 * (u₂ - 2u₃ + u₄)^2 + 1 / 4 * (u₂ - u₄)^2
    β₃ = 13 / 12 * (u₃ - 2u₄ + u₅)^2 + 1 / 4 * (3u₃ - 4u₄ + u₅)^2
    w̃₁ = (1 / 10) / (β₁ + ϵ)^2
    w̃₂ = (3 / 5) / (β₂ + ϵ)^2
    w̃₃ = (3 / 10) / (β₃ + ϵ)^2
    s = w̃₁ + w̃₂ + w̃₃
    w₁, w₂, w₃ = w̃₁ / s, w̃₂ / s, w̃₃ / s
    w₁ * v₁ + w₂ * v₂ + w₃ * v₃
end
