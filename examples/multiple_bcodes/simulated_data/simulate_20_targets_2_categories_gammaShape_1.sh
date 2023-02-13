java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed 1234 simulate_tree.xml

for seed in `seq 1 10`

do

java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed $seed -D "clock=0.0343711018512661" simulate_2_categories_small.xml

done

for seed in `seq 1 10`

do

java -jar ../../../../out/artifacts/beast_gestalt_jar/beast_gestalt.jar -seed $seed -D "clock=0.1656288981487339" simulate_2_categories_small.xml

done


