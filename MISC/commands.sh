#run for the 10k dataset gene tree 1
# for n in {50,100,200,500,1000,10000}; do
n="1000"
k=5
for rep in 0{1..9} 10; do 
	for s in {1..10}; do
	sbatch run_opt.sh $n $rep $s $k 0
done
done
# done

n="100"
k=5
for rep in 0{1..9} 10; do 
	for s in {1..10}; do
	sbatch --exclude=n0005,n0020 run_opt.sh $n $rep $s $k 0
done
done

python find_sol.py -i ../quartet_SPR/Dataset/10k_dataset/100/01/gene_tree_1 -s 10 -k 5 -n 1 -l 10000 -o ./Datasets/10k/100/01/gene_tree_1_k5_s10_noise1_l10000

grep "Jaccard" Datasets/10k/200/*/gene_tree_1_k5_s*_noise0/log.txt

grep "Jaccard" Datasets/10k/100/*/gene_tree_1_k5_s*_noise1_l10000/log.txt

rm ./Datasets/10k/runtime_k5.txt
for n in {50,100,200,500,1000}; do
for rep in 0{1..9} 10; do 
for s in {1..10}; do
	python compute_runtime.py -i Datasets/10k/$n/$rep/gene_tree_1_k5_s${s}_noise0 | awk -v var1="$n" -v var2="$rep" -v var3="$s" -v var4="5" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' >> ./Datasets/10k/runtime_k5.txt
done
done
done

rm ./Datasets/10k/runtime_n100.txt
# for n in {50,100,200,500,1000}; do
for rep in 0{1..9} 10; do 
for s in {1..10}; do
for k in {2..10}; do
	python compute_runtime.py -i Datasets/10k/100/$rep/gene_tree_1_k${k}_s${s}_noise0 | awk -v var1="100" -v var2="$rep" -v var3="$s" -v var4="$k" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' >> ./Datasets/10k/runtime_n100.txt
done
done
done

rm ./Datasets/10k/pr_n100.txt
n=100
for rep in 0{1..9} 10; do 
for s in {1..10}; do
	python compute_runtime.py -i Datasets/10k/$n/$rep/gene_tree_1_k5_s${s}_f1_noise0 | awk -v var1="$n" -v var2="$rep" -v var3="$s" -v var4="5" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' >> ./Datasets/10k/pr_n100.txt
done
done


rm ./Datasets/10k/jaccard_n100.txt
n=100
for k in {2..10}; do
for rep in 0{1..9} 10; do 
for s in {1..10}; do
	grep "Jaccard" Datasets/10k/gene_tree_1_${n}_${rep}_k${k}_s${s}.log | awk -v var1="$n" -v var2="$rep" -v var3="$s" -v var4="$k" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $2}' >> ./Datasets/10k/jaccard_n100.txt
done
done
done

rm ./Datasets/10k/runtime_n100.txt
n=100
for k in {2..10}; do
for rep in 0{1..9} 10; do 
for s in {1..10}; do
	grep "Optimization" Datasets/10k/gene_tree_1_${n}_${rep}_k${k}_s${s}.log | awk -v var1="$n" -v var2="$rep" -v var3="$s" -v var4="$k" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $3}' >> ./Datasets/10k/runtime_n100.txt
done
done
done

#generate krepp input
# for t in {1..7}; do
for t in {1068,202,1322}; do
for k in {3,5,10}; do
for s in {1..5}; do
	rm -r 10k_dataset/tree_${t}_k${k}_s${s}_uniform
	rm -r 10k_dataset/tree_${t}_k${k}_s${s}_exp
	mkdir -p 10k_dataset/tree_${t}_k${k}_s${s}_uniform
	mkdir -p 10k_dataset/tree_${t}_k${k}_s${s}_exp
	python generate_krepp_input.py -i 10k_dataset/tree_$t.tre -s $s -k $k -o 10k_dataset/tree_${t}_k${k}_s${s}_uniform -p "uniform"
	python generate_krepp_input.py -i 10k_dataset/tree_$t.tre -s $s -k $k -o 10k_dataset/tree_${t}_k${k}_s${s}_exp -p "exp"
done
done
done

python generate_krepp_input.py -i 10k_dataset/tree_7.tre -s 1 -k 5 -o 10k_dataset -p "uniform"

#run for 10k dataset

# for t in {4,7} {8..13}; do
for t in {4,7,15,17,1322}; do
# nw_reroot ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre > ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre.rerooted
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
sbatch --exclude=n0020,n0005 run_opt_10k.sh $t $k $s $u
done
done
done
done

for t in {1068,202,1322}; do
# nw_reroot ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre > ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre.rerooted
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
# wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_distances.txt
cp ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_distances.txt ../Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u
done
done
done
done


t=4
k=10
s=1
# u="exp"
for u in {uniform,exp}; do
sbatch --exclude=n0020,n0005 run_opt_10k.sh $t $k $s $u
done

# for t in {15,17,1322}; do
# t=1322
for t in {4,7,15,17,1322}; do
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
# wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist/metrics.txt
# awk '{print $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/metrics.txt
# wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_dist/metrics.txt
# wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_distances_pr-h14_dth4.txt 
# wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_mash_metrics.txt
# awk '{print $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/ani_blen_tree_estimated_h14_dth4_dist/metrics.txt
wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/ani_blen_tree_estimated_corrected_h14_dth4_dist/metrics.txt
done
done
done
done

t=15
k=10
for s in {1..5}; do
for u in {uniform,exp}; do
	grep "Optimization Runtime" ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/log.txt
done
done

# rm ./Demixing_Data/10k_dataset/all_metrics_true.txt
# rm ./Demixing_Data/10k_dataset/all_metrics_true_mash.txt
# rm ./Demixing_Data/10k_dataset/all_metrics_est.txt
# rm ./Demixing_Data/10k_dataset/all_metrics_est_mash.txt
rm ./Demixing_Data/10k_dataset/all_metrics_est_mash_corrected.txt
for t in {4,7,15,17,1322}; do
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
# awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist/metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_est.txt
# awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/ani_blen_tree_estimated_h14_dth4_dist/metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_est_mash.txt
awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/ani_blen_tree_estimated_corrected_h14_dth4_dist/metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_est_mash_corrected.txt
# awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_true.txt
# awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_dist/metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_true_mash.txt
done
done
done
done

t=1322
k=5
s=1
u="exp"
python find_sol_10k.py -i ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre.rerooted -d ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u -e 1 -o ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_mash_dist/log.txt
python find_sol_10k.py -i ./Demixing_Data/10k_dataset/tree_${t}-ani_blen.tre.rerooted -d ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u -e 0 -o ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/log.txt

rm -r ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist
mkdir -p ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist
python find_sol_10k.py -i ./Demixing_Data/10k_dataset/tree_$t.tre -d ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u -e 0 -o ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/log.txt
python compute_metrics.py -i ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/metrics.txt

rm -r ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist
mkdir -p ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist
python find_sol_10k.py -i ./Demixing_Data/10k_dataset/tree_$t.tre -d ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u -e 1 -o ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist/log.txt
python compute_metrics.py -i ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist > ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist/metrics.txt


for t in ; do
t=1322
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
sbatch run_metric_10k.sh $t $k $s $u
done
done
done
done

for t in {4,7} {14..17}; do
t=1322
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_metrics.txt
wc -l ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_mash_metrics.txt
done
done
done
done


rm ./Demixing_Data/10k_dataset/all_metrics_voltka.txt
rm ./Demixing_Data/10k_dataset/all_metrics_voltka_mash.txt
for t in {4,7,15,17,1322}; do
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_voltka.txt
awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_mash_metrics.txt >> ./Demixing_Data/10k_dataset/all_metrics_voltka_mash.txt
done
done
done
done

t=4
k=3
s=1
u="exp"
python compute_metrics_voltka.py -t Demixing_Data/10k_dataset/tree_$t.tre -i Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u > Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/voltka_metrics.txt


#get metrics
rm ./Demixing_Data/10k_dataset/true_dist_res.txt
for t in {1,4,7}; do
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/true_dist/metrics.txt >> ./Demixing_Data/10k_dataset/true_dist_res.txt
done
done
done
done

rm ./Demixing_Data/10k_dataset/estimated_h14_dth4_dist_res.txt
for t in {1,4,7}; do
for k in {3,5,10}; do
for s in {1..5}; do
for u in {uniform,exp}; do
awk -v var1="$t" -v var2="$k" -v var3="$s" -v var4="$u" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' ./Demixing_Data/10k_dataset/tree_${t}_k${k}_s${s}_$u/estimated_h14_dth4_dist/metrics.txt~ >> ./Demixing_Data/10k_dataset/estimated_h14_dth4_dist_res.txt
done
done
done
done
#Bio Data

#bees-bosseret -> done
d="bees-bosseret"
t="bees-bosseret.tre"

#mammals-song -> done
d="mammals-song"
t="mammals-song-castlespro.tre"

#birds-jarvis -> done
d="birds-jarvis"
t="birds-jarvis.tre.rerooted"

#tilapia-Ciezarek -> done
d="tilapia-Ciezarek"
t="tilapia-Ciezarek.tre.rerooted"

#1kp -> done
d="1kp"
t="1kp-concat-fig2.tre"

#pancrustacean-Bernot -> done
d="pancrustacean-Bernot"
t="pancrustacean-Bernot.tre.rerooted"

#beetles-Johnson -> done
d="beetles-Johnson"
t="beetles-Johnson.tre"

#fish-Troyer -> done
d="fish-Troyer"
t="fish-Troyer.tre.rerooted"

#hemipteroid-johnson -> done
d="hemipteroid-johnson"
t="hemipteroid-johnson.tre.rerooted"

#mammals-foley -> done
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

#fish_hughes -> done
d="fish_hughes"
t="fish_hughes_1105examl.tre.rerooted"

#birds-stiller -> done
d="birds-stiller"
t="birds-stiller-castlepro.tre"

for k in {2,3,5,7,10}; do
for s in {1..10}; do
for n in {0,1,2}; do
sbatch --exclude=n0020 run_opt_bio.sh Demixing_Data/biotrees/$d $t $s $k $n
done
done
done


grep "Jaccard" ./Demixing_Data/biotrees/$d/${t}_k*_s*_noise0/log.txt | wc -l

python find_sol.py -i ./Demixing_Data/biotrees/$d/${t} -s 3 -k 3 -n 0 -o ./Demixing_Data


#unifrac results
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

sbatch run_metrics_bio.sh $d $t

rm ./Demixing_Data/biotrees/all_unifrac_res.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot","beetles-Johnson","fish-Troyer","hemipteroid-johnson","mammals-foley","fish_hughes","birds-stiller"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees/$d/unifrac_res.txt >> ./Demixing_Data/biotrees/all_unifrac_res.txt
done

rm ./Demixing_Data/biotrees/all_runtime_res.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot","beetles-Johnson","fish-Troyer","hemipteroid-johnson","mammals-foley","fish_hughes","birds-stiller"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees/$d/runtime_res.txt >> ./Demixing_Data/biotrees/all_runtime_res.txt
done


rm ./Demixing_Data/biotrees/$d/unifrac_res_new.txt
for k in {3,5,7,10}; do
# for k in {2,3,5,7,10}; do
# for s in {1..10}; do
for s in {1..5}; do
paste  <(python compute_metrics.py -i ./Demixing_Data/biotrees/$d/${t}_k${k}_s${s}_noise0) <(python compute_metrics.py -i ./Demixing_Data/biotrees/$d/${t}_k${k}_s${s}_noise0_new) | awk -v var1="$k" -v var2="$s" -v var3="0" -v var4="0" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' >> ./Demixing_Data/biotrees/$d/unifrac_res_new.txt
done
done

for k in {3,5,7,10}; do
# for k in {2,3,5,7,10}; do
# for s in {1..10}; do
for s in {1..5}; do
for l in {100,1000,10000}; do
paste <(python compute_metrics.py -i ./Demixing_Data/biotrees/$d/${t}_k${k}_s${s}_noise1_l${l})  <(python compute_metrics.py -i ./Demixing_Data/biotrees/$d/${t}_k${k}_s${s}_noise1_l${l}_new) | awk -v var1="$k" -v var2="$s" -v var3="1" -v var4="$l" '{print var1 "\t" var2 "\t" var3 "\t" var4 "\t" $0}' >> ./Demixing_Data/biotrees/$d/unifrac_res_new.txt
done
done
done

d="mammals-foley"
mkdir ./Demixing_Data/biotrees/$d
scp sarasti@calab-fe.ucsd.edu:../../../../scratch00/sarasti/demixing/Demixing_Data/biotrees/$d/unifrac_res.txt ./Demixing_Data/biotrees/$d
scp sarasti@calab-fe.ucsd.edu:../../../../scratch00/sarasti/demixing/Demixing_Data/biotrees/$d/jaccard_res.txt ./Demixing_Data/biotrees/$d


#frogs-feng -> check 1 and 2
for k in {7,10}; do
for s in {1..10}; do
sbatch run_opt_bio.sh Demixing_Data/biotrees/frogs-feng frogs-feng.tre.rerooted $s $k 1
done
done

grep "Jaccard" ./Demixing_Data/biotrees/frogs-feng/frogs-feng.tre.rerooted_k*_s*_noise0/log.txt


python find_sol.py -i ./Demixing_Data/biotrees/frogs-feng/frogs-feng.tre.rerooted -s 4 -k 7 -n 1 -l 100 -o ./Demixing_Data

#turtle-thomas -> done 
for k in {2,3,5,7,10}; do
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/turtle-thomas turtle-thomas.tre.rerooted $s $k 0
done
done

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/turtle-thomas turtle-thomas.tre.rerooted $s 0
done

grep "Jaccard" ./Demixing_Data/biotrees/turtle-thomas/turtle-thomas.tre.rerooted_k*_s*_noise0/log.txt

#plants-onekp -> done
for k in {2,3,5,7,10}; do
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/plants-onekp plants-onekp-castlespro.tre $s $k 0
done
done

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/plants-onekp plants-onekp-castlespro.tre $s 0
done

grep "Jaccard" ./Demixing_Data/biotrees/plants-onekp/plants-onekp-castlespro.tre_k*_s*_noise0/log.txt

python find_sol.py -i ./Demixing_Data/biotrees/plants-onekp/plants-onekp-castlespro.tre -s 1 -k 10 -n 0 -o ./Demixing_Data

#fungi-li

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/fungi-li fungi-li.tre.rerooted $s 0
done

#subosciences-harey

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/subosciences-harey subosciences-harey.tre.rerooted $s 0
done

#mammals-Upham

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/mammals-Upham mammals-Upham.tre.rerooted $s 0
done

#bees-Piskulich
for k in {2,3,5,7,10}; do
for s in {1..10}; do
sbatch run_opt_bio.sh Demixing_Data/biotrees/bees-Piskulich bees-Piskulich.tre.rerooted $s $k 0
done
done

#expanse
for s in {1..10}; do
sbatch run_opt_bio.sh ./Demixing_Data/biotrees/bees-Piskulich bees-Piskulich.tre.rerooted $s 0
done

grep "Jaccard" ./Demixing_Data/biotrees/bees-Piskulich/bees-Piskulich.tre.rerooted_k*_s*_noise0/log.txt

python find_sol.py -i ./Demixing_Data/biotrees/bees-Piskulich/bees-Piskulich.tre.rerooted -s 1 -k 2 -n 0 -o ./Demixing_Data


##Compare with baseline

#fixed k hill climbing

#bees-bosseret -> done
d="bees-bosseret"
t="bees-bosseret.tre"

#mammals-song -> done
d="mammals-song"
t="mammals-song-castlespro.tre"

#birds-jarvis -> done
d="birds-jarvis"
t="birds-jarvis.tre.rerooted"

#tilapia-Ciezarek -> check
d="tilapia-Ciezarek"
t="tilapia-Ciezarek.tre.rerooted"

#1kp -> done
d="1kp"
t="1kp-concat-fig2.tre"

#pancrustacean-Bernot -> done
d="pancrustacean-Bernot"
t="pancrustacean-Bernot.tre.rerooted"

#beetles-Johnson -> done
d="beetles-Johnson"
t="beetles-Johnson.tre"

#fish-Troyer -> run mter
d="fish-Troyer"
t="fish-Troyer.tre.rerooted"

#hemipteroid-johnson -> run met
d="hemipteroid-johnson"
t="hemipteroid-johnson.tre.rerooted"

#mammals-foley -> run mter
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

#fish_hughes -> run mter
d="fish_hughes"
t="fish_hughes_1105examl.tre.rerooted"

#frogs-feng -> done
d="frogs-feng"
t="frogs-feng.tre.rerooted"

#birds-stiller -> done
d="birds-stiller"
t="birds-stiller-castlepro.tre"

for n in {0,1}; do
# for k in {2,3,5,7,10}; do
# for s in {1..10}; do
s=5
n=0
for k in {7,10}; do
sbatch --exclude=n0020 run_comp_bio.sh Demixing_Data/biotrees/$d/$t $s $k $n Demixing_Data/biotrees_compare_methods/$d/${t}_k${k}_s${s}_f1_noise$n "hill"
done
done
done

grep "Jaccard" ./Demixing_Data/biotrees_compare_methods/$d/${t}_k*_s*_f1_noise0/log.txt

#compute metrics
sbatch run_metrics_comp_bio.sh $d $t

rm ./Demixing_Data/biotrees/all_unifrac_res_f1.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot","beetles-Johnson","fish-Troyer","hemipteroid-johnson","mammals-foley","fish_hughes","birds-stiller"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees_compare_methods/$d/unifrac_res_f1.txt >> ./Demixing_Data/biotrees/all_unifrac_res_f1.txt
done

#fixed-k exhaustive
#bees-bosseret -> done
d="bees-bosseret"
t="bees-bosseret.tre"

#mammals-song -> done
d="mammals-song"
t="mammals-song-castlespro.tre"

#birds-jarvis -> done
d="birds-jarvis"
t="birds-jarvis.tre.rerooted"

#tilapia-Ciezarek -> done
d="tilapia-Ciezarek"
t="tilapia-Ciezarek.tre.rerooted"

#1kp -> done
d="1kp"
t="1kp-concat-fig2.tre"

#pancrustacean-Bernot -> done
d="pancrustacean-Bernot"
t="pancrustacean-Bernot.tre.rerooted"

#beetles-Johnson -> up to here
d="beetles-Johnson"
t="beetles-Johnson.tre"

#fish-Troyer
d="fish-Troyer"
t="fish-Troyer.tre.rerooted"

#hemipteroid-johnson
d="hemipteroid-johnson"
t="hemipteroid-johnson.tre.rerooted"

#mammals-foley
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

#fish_hughes
d="fish_hughes"
t="fish_hughes_1105examl.tre.rerooted"

for n in {0,1,2}; do
for k in {2,3}; do
for s in {1..10}; do
k=2
s=10
n=1
sbatch run_comp_bio.sh Demixing_Data/biotrees/$d/$t $s $k $n Demixing_Data/biotrees_compare_methods/$d/${t}_k${k}_s${s}_f1_exh_noise$n "exhaustive"
done
done
done

grep "Jaccard" ./Demixing_Data/biotrees_compare_methods/$d/${t}_k*_s*_f1_exh_noise0/log.txt

#compute metrics
sbatch run_metrics_comp_bio_closest.sh $d $t "exh"

rm ./Demixing_Data/biotrees/all_unifrac_res_f1_exh.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees_compare_methods/$d/unifrac_res_f1_exh.txt >> ./Demixing_Data/biotrees/all_unifrac_res_f1_exh.txt
done

#fixed-k closest
#bees-bosseret -> done
d="bees-bosseret"
t="bees-bosseret.tre"

#mammals-song -> done
d="mammals-song"
t="mammals-song-castlespro.tre"

#birds-jarvis -> done
d="birds-jarvis"
t="birds-jarvis.tre.rerooted"

#tilapia-Ciezarek -> done
d="tilapia-Ciezarek"
t="tilapia-Ciezarek.tre.rerooted"

#1kp -> done
d="1kp"
t="1kp-concat-fig2.tre"

#pancrustacean-Bernot -> done
d="pancrustacean-Bernot"
t="pancrustacean-Bernot.tre.rerooted"

#beetles-Johnson -> done
d="beetles-Johnson"
t="beetles-Johnson.tre"

#fish-Troyer -> done
d="fish-Troyer"
t="fish-Troyer.tre.rerooted"

#hemipteroid-johnson -> done
d="hemipteroid-johnson"
t="hemipteroid-johnson.tre.rerooted"

#mammals-foley -> done
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

#fish_hughes -> done
d="fish_hughes"
t="fish_hughes_1105examl.tre.rerooted"

#frogs-feng -> done
d="frogs-feng"
t="frogs-feng.tre.rerooted"

#birds-stiller -> done
d="birds-stiller"
t="birds-stiller-castlepro.tre"

# for n in {0,1,2}; do
for k in {2,3,5,7,10}; do
for s in {1..10}; do
n=2
k=10
sbatch run_comp_bio_closest.sh Demixing_Data/biotrees/$d/$t $s $k $n Demixing_Data/biotrees_compare_methods/$d/${t}_k${k}_s${s}_f1_clos_noise$n
done
done
done

grep "Jaccard" ./Demixing_Data/biotrees_compare_methods/$d/${t}_k*_s*_f1_clos_noise0/log.txt

#compute metrics
sbatch run_metrics_comp_bio_closest.sh $d $t "clos"

rm ./Demixing_Data/biotrees/all_unifrac_res_f1_clos.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot","beetles-Johnson","fish-Troyer","hemipteroid-johnson","mammals-foley","fish_hughes","birds-stiller"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees_compare_methods/$d/unifrac_res_f1_clos.txt >> ./Demixing_Data/biotrees/all_unifrac_res_f1_clos.txt
done

#fixed-k closest iterative

#bees-bosseret -> done
d="bees-bosseret"
t="bees-bosseret.tre"

#mammals-song -> done
d="mammals-song"
t="mammals-song-castlespro.tre"

#birds-jarvis -> done
d="birds-jarvis"
t="birds-jarvis.tre.rerooted"

#tilapia-Ciezarek -> done
d="tilapia-Ciezarek"
t="tilapia-Ciezarek.tre.rerooted"

#1kp -> done
d="1kp"
t="1kp-concat-fig2.tre"

#pancrustacean-Bernot -> done
d="pancrustacean-Bernot"
t="pancrustacean-Bernot.tre.rerooted"

#beetles-Johnson -> done
d="beetles-Johnson"
t="beetles-Johnson.tre"

#fish-Troyer -> done
d="fish-Troyer"
t="fish-Troyer.tre.rerooted"

#hemipteroid-johnson -> done
d="hemipteroid-johnson"
t="hemipteroid-johnson.tre.rerooted"

#mammals-foley -> done
d="mammals-foley"
t="mammals-foley-concat-heutral-241.tre"

#fish_hughes -> done
d="fish_hughes"
t="fish_hughes_1105examl.tre.rerooted"

#frogs-feng -> done
d="frogs-feng"
t="frogs-feng.tre.rerooted"

#birds-stiller -> done
d="birds-stiller"
t="birds-stiller-castlepro.tre"

for n in {0,1,2}; do
for k in {2,3,5,7,10}; do
for s in {1..10}; do
# for k in {7,10}; do
# n=2
# k=10
sbatch run_comp_bio_closest.sh Demixing_Data/biotrees/$d/$t $s $k $n Demixing_Data/biotrees_compare_methods/$d/${t}_k${k}_s${s}_f1_clos_iter_noise$n
done
done
done

grep "Jaccard" ./Demixing_Data/biotrees_compare_methods/$d/${t}_k*_s*_f1_clos_iter_noise0/log.txt

#compute metrics
sbatch run_metrics_comp_bio_closest.sh $d $t "clos_iter"

rm ./Demixing_Data/biotrees/all_unifrac_res_f1_clos_iter.txt
for d in {"bees-bosseret","mammals-song","birds-jarvis","tilapia-Ciezarek","1kp","pancrustacean-Bernot","beetles-Johnson","fish-Troyer","hemipteroid-johnson","mammals-foley","fish_hughes","birds-stiller"}; do
	awk -v var1="$d" '{print var1 "\t" $0}' ./Demixing_Data/biotrees_compare_methods/$d/unifrac_res_f1_clos_iter.txt >> ./Demixing_Data/biotrees/all_unifrac_res_f1_clos_iter.txt
done

#TreeShrink Dataset
for t in {"Plants","Mammals","Insects","Frogs","XenCannon","XenRouse"}; do
	mkdir -p ./Datasets/$t
	head -n1 ../tree_matching/Datasets/TreeShrink_data/$t/unfiltered.trees.rerooted > ./Datasets/$t/unfiltered.trees.rerooted_1
done

k=5
for t in {"Plants","Mammals","Insects","Frogs","XenCannon","XenRouse"}; do
	for s in {1..10}; do
	for l in {100,1000,10000}; do
	sbatch run_opt_bio.sh $t $s $k 2 $l
	# sbatch run_opt_bio.sh $t $s $k 0
	done
done
done

for t in {"Plants","Mammals","Insects","Frogs","XenCannon","XenRouse"}; do
grep "Jaccard" Datasets/Mammals/unfiltered_trees_rerooted_1_k5_s*_noise0/log.txt | wc -l
done

for l in {100,1000,10000}; do
for t in {"Plants","Mammals","Insects","Frogs","XenCannon","XenRouse"}; do
grep "Jaccard" Datasets/$t/unfiltered_trees_rerooted_1_k5_s*_noise1_l$l/log.txt
done
done

mkdir -p ./Datasets/$1/unfiltered_trees_rerooted_1_k$3_s$2_noise$4
python find_sol.py -i ./Datasets/$1/unfiltered.trees.rerooted_1 -s $2 -k $3 -n $4 -o ./Datasets/$1/unfiltered_trees_rerooted_1_k$3_s$2_noise$4 > ./Datasets/$1/unfiltered_trees_rerooted_1_k$3_s$2_noise$4/log.txt

python find_sol.py -i ./Datasets/Mammals/unfiltered.trees.rerooted_1 -s 1 -k 5 -n 0 -o ./Datasets/Mammals/unfiltered_trees_rerooted_1_k5_s1_noise0


