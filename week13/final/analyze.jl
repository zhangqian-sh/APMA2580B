using Statistics
using DelimitedFiles
using JLD2
using PrettyTables, Printf

setup = ARGS[1]
t = parse(Float64, ARGS[2])
# exclude_shock = parse(Bool, ARGS[3])
exclude_shock = true
path = "./misc/output/$(setup)-test/t=$(t)/"

l1(x, y) = mean(abs.(x - y))
l2(x, y) = sqrt(mean((x - y) .^ 2))
l∞(x, y) = maximum(abs.(x - y))

# N = [20, 40, 80, 160, 320, 640]
N = [20, 40, 80, 160, 320, 640]
data = zeros(Float64, (length(N), 6))
const newaxis = [CartesianIndex()]
# compute error
if setup == "density_wave" || setup == "burgers_like"
    for i in eachindex(N[1:end])
        # u, u_exact = readdlm(path * "u_N=$(N[i]).csv", ','), readdlm(path * "u_exact_N=$(N[i]).csv", ',')
        u = jldopen(path * "u_N=$(N[i]).jld2", "r")
        u_exact = jldopen(path * "u_exact_N=$(N[i]).jld2", "r")
        if setup == "burgers_like"
            u = vcat(u["rho"][newaxis, :, :], u["v"][newaxis, :, :], u["w"][newaxis, :, :])
            u_exact = vcat(u_exact["rho"][newaxis, :, :], u_exact["v"][newaxis, :, :], u_exact["w"][newaxis, :, :])
            if exclude_shock
                x = y = LinRange(0, 2√2π, N[i])
                for idx_i = 1:N[i]
                    for idx_j = 1:N[i]
                        ξ = (x[idx_j] + y[idx_i]) / √2
                        if π < ξ < π + 0.5π || 3π < ξ < 3π + 0.5π
                            u[:, idx_i, idx_j] .= NaN
                        end
                    end
                end
            end
        else
            u = vcat(u["rho"][newaxis, :, :], u["v"][newaxis, :, :], u["w"][newaxis, :, :], u["p"][newaxis, :, :])
            u_exact = vcat(u_exact["rho"][newaxis, :, :], u_exact["v"][newaxis, :, :], u_exact["w"][newaxis, :, :], u_exact["p"][newaxis, :, :])
        end
        nanmean(x) = mean(filter(!isnan, x))
        nanmanximum(x) = maximum(filter(!isnan, x))
        nan_l1(x, y) = nanmean(abs.(x - y))
        nan_l2(x, y) = sqrt(nanmean((x - y) .^ 2))
        nan_l∞(x, y) = nanmanximum(abs.(x - y))
        data[i, 1], data[i, 3], data[i, 5] = nan_l1(u, u_exact), nan_l2(u, u_exact), nan_l∞(u, u_exact)
    end
elseif setup == "Lax" || setup == "Sod"
    for i in eachindex(N[1:end-1])
        u, u_exact = readdlm(path * "u_N=$(N[i]).csv", ','), readdlm(path * "u_N=$(N[i+1]).csv", ',')
        u_exact = u_exact[:, 1:2:end]
        data[i, 1], data[i, 3], data[i, 5] = l1(u, u_exact), l2(u, u_exact), l∞(u, u_exact)
    end
end
# compute order
for i = 2:length(N)
    data[i, 2] = log2(data[i-1, 1] / data[i, 1])
    data[i, 4] = log2(data[i-1, 3] / data[i, 3])
    data[i, 6] = log2(data[i-1, 5] / data[i, 5])
end

formatters = (v, i, j) -> mod(j, 2) == 1 ? ft_printf("%.2e")(v, i, j) : ft_printf("%.2f")(v, i, j)
pretty_table(data, header=["L_1"; "Order"; "L_2"; "Order"; "L_∞"; "Order"], formatters=(formatters, ft_nonothing))
path = "misc/table/$(setup)-test/"
mkpath(path)
suffix = exclude_shock ? "_es" : ""
writedlm(path * "t=$(t)" * suffix * ".csv", data, ",")

# turn data into tex
# println(data)
# println(size(data))
lines = Vector{String}()
for row_idx in axes(data, 1)
    l = data[row_idx, :]
    nums = Vector{String}()
    push!(nums, "$(2^(row_idx)*10)")
    # println(row_idx, " ", l)
    for (col_idx, num) in enumerate(l)
        # println(col_idx, " ", num)
        if mod(col_idx, 2) == 1
            s = @sprintf("%.2e", num)
            sgn = log10(num) < 0 ? "-" : ""
            s = "\\($(s[1:4])\\times10^{$(sgn)$(parse(Int, s[7:end]))}\\)"
        else
            s = row_idx > 1 ? @sprintf("%.2f", num) : "-"
        end
        push!(nums, s)
        # println(s, " ", nums)
    end
    line = join(nums, " & ")
    # println(line)
    push!(lines, line)
end

lines = join(lines, "\\\\\n") * "\\\\\n"

open(path * "t=$(t)" * suffix * ".tex", "w") do io
    write(io, lines)
end;