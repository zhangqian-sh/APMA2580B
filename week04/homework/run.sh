flux=LaxWendroff
t=0.3
for N_x in 20 40 80 160 320 640
do
    julia main.jl $t $N_x
done
julia analyze.jl $flux $t false> ./misc/table/"${flux}_t=${t}.txt"
t=2.0
for N_x in 20 40 80 160 320 640
do
    julia main.jl $t $N_x
done
julia analyze.jl $flux $t false> ./misc/table/"${flux}_t=${t}.txt"
julia analyze.jl $flux $t true> ./misc/table/"${flux}_t=${t}_es.txt"
