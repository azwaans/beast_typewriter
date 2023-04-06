#this script is used to simulate 20 targetBCs with a clock rate from a 2 categories site model with GammaShape = 1 (10 each).
#1st, simulate a tree
java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed 1234 simulate_tree.xml

#2nd, simulate 10 targetBC alignments for each of the 2 category clock values
for seed in `seq 1 10`

do

java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed $seed -D "clock=0.0343711018512661" simulate_alignment_fixed_tree_clock_input.xml

done

for seed in `seq 1 10`

do

java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed $seed -D "clock=0.1656288981487339" simulate_alignment_fixed_tree_clock_input.xml

done


