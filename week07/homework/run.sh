limiter=minimod_limiter
t=0.3
for N_x in 20 40 80 160 320 640
do
    julia main.jl $t $N_x
done
julia analyze.jl $limiter $t false> ./misc/table/"${limiter}_t=${t}.txt"
t=2.0
for N_x in 20 40 80 160 320 640
do
    julia main.jl $t $N_x
done
julia analyze.jl $limiter $t false> ./misc/table/"${limiter}_t=${t}.txt"
julia analyze.jl $limiter $t true> ./misc/table/"${limiter}_t=${t}_es.txt"
