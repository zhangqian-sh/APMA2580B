using LinearAlgebra

function WENO5_reconstruction(u₁, u₂, u₃, u₄, u₅)::Float64
    ϵ = 1e-6
    v₁ = 1 / 3 * u₁ - 7 / 6 * u₂ + 11 / 6 * u₃
    v₂ = -1 / 6 * u₂ + 5 / 6 * u₃ + 1 / 3 * u₄
    v₃ = 1 / 3 * u₃ + 5 / 6 * u₄ - 1 / 6 * u₅
    β₁ = 13 / 12 * (u₁ - 2 * u₂ + u₃)^2 + 1 / 4 * (u₁ - 4 * u₂ + 3 * u₃)^2
    β₂ = 13 / 12 * (u₂ - 2 * u₃ + u₄)^2 + 1 / 4 * (u₂ - u₄)^2
    β₃ = 13 / 12 * (u₃ - 2 * u₄ + u₅)^2 + 1 / 4 * (3 * u₃ - 4 * u₄ + u₅)^2
    # β₁, β₂, β₃ = 1, 1, 1
    w̃₁ = (1 / 10) / (β₁ + ϵ)^2
    w̃₂ = (3 / 5) / (β₂ + ϵ)^2
    w̃₃ = (3 / 10) / (β₃ + ϵ)^2
    s = w̃₁ + w̃₂ + w̃₃
    w₁, w₂, w₃ = w̃₁ / s, w̃₂ / s, w̃₃ / s
    w₁ * v₁ + w₂ * v₂ + w₃ * v₃
end