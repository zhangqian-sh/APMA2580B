t=0.3
for M in 0.0 1.0 5.0 10.0
do 
    for N_x in 20 40 80 160 320 640
    do
        julia main.jl $t $N_x $M
    done
    mkdir -p "./misc/table/M=${M}/"
    julia analyze.jl $M $t false> ./misc/table/"M=${M}/t=${t}.txt"
done

t=2.0
for M in 0.0 1.0 5.0 10.0
do 
    for N_x in 20 40 80 160 320 640
    do
        julia main.jl $t $N_x $M
    done
    mkdir -p "./misc/table/M=${M}/"
    julia analyze.jl $M $t false> ./misc/table/"M=${M}/t=${t}.txt"
    julia analyze.jl $M $t true> ./misc/table/"M=${M}/t=${t}_es.txt"
done
