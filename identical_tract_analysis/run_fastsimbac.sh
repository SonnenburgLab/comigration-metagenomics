num_samples=100
# each sub pop is half of the total pop size
sub_pop_size=50
seq_length=1000000
recomb_rate=0.005
mut_rate=0.01
recomb_length=10000
save_path="fastsimbac_res/l=$recomb_length"

# program path
fastsimbac_path="/Users/Device6/Documents/Research/bgoodlab/bioinformatics/fastsimbac/fastSimBac_mac/fastsimbac"

# sim reps
# for t in 0 0.01 0.03 0.05 0.07 0.1 0.3 0.5 0.7 1 3
for t in 0.001 0.003 0.005 0.007
do
    time $fastsimbac_path $num_samples $seq_length -r $recomb_rate $recomb_length -I 2 $sub_pop_size $sub_pop_size 0 -t $mut_rate -ej $t 2 1 1>$save_path/t=$t.txt
    echo "t=$t done"
done
