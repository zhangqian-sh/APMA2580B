using Statistics
using DelimitedFiles
using PrettyTables, Printf

M, t = parse(Float64, ARGS[1]), parse(Float64, ARGS[2])
exclude_shock = parse(Bool, ARGS[3])
path = "./misc/output/M=$(M)/t=$(t)/"

l1(x, y) = mean(abs.(x - y))
l2(x, y) = sqrt(mean((x - y) .^ 2))
l∞(x, y) = maximum(abs.(x - y))

N = [20, 40, 80, 160, 320, 640]
data = zeros(Float64, (length(N), 6))
# compute error
for i in eachindex(N)
    u, u_exact = readdlm(path * "u_N=$(N[i]).csv"), readdlm(path * "u_exact_N=$(N[i]).csv")
    if exclude_shock
        x_interface = LinRange(0, 2π, N[i] + 1)
        x = x_interface[1:end-1] .+ π / N[i]
        shock_position = π + 2 / 3
        u, u_exact = u[abs.(x .- shock_position).>=1], u_exact[abs.(x .- shock_position).>=1]
    end
    data[i, 1], data[i, 3], data[i, 5] = l1(u, u_exact), l2(u, u_exact), l∞(u, u_exact)
end
# compute order
for i = 2:length(N)
    data[i, 2] = log2(data[i-1, 1] / data[i, 1])
    data[i, 4] = log2(data[i-1, 3] / data[i, 3])
    data[i, 6] = log2(data[i-1, 5] / data[i, 5])
end

formatters = (v, i, j) -> mod(j, 2) == 1 ? ft_printf("%.2e")(v, i, j) : ft_printf("%.2f")(v, i, j)
pretty_table(data, header=["L_1"; "Order"; "L_2"; "Order"; "L_∞"; "Order"], formatters=(formatters, ft_nonothing))
path = "misc/table/M=$(M)/"
mkpath(path)
suffix = exclude_shock ? "_es" : ""
writedlm(path * "t=$(t)" * suffix * ".csv", data, ",")

# turn data into tex
println(data)
println(size(data))
lines = Vector{String}()
for row_idx in axes(data, 1)
    l = data[row_idx, :]
    nums = Vector{String}()
    push!(nums, "$(2^(row_idx)*10)")
    println(row_idx, " ", l)
    for (col_idx, num) in enumerate(l)
        println(col_idx, " ", num)
        if mod(col_idx, 2) == 1
            s = @sprintf("%.2e", num)
            sgn = log10(num) < 0 ? "-" : ""
            s = "\\($(s[1:4])\\times10^{$(sgn)$(s[end])}\\)"
        else
            s = row_idx > 1 ? @sprintf("%.2f", num) : "-"
        end
        push!(nums, s)
        println(s, " ", nums)
    end
    line = join(nums, " & ")
    println(line)
    push!(lines, line)
end

lines = join(lines, "\\\\\n") * "\\\\\n"

open(path * "t=$(t)" * suffix * ".tex", "w") do io
    write(io, lines)
end;