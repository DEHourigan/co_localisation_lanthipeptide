# this script was used on the segments that were hits from hmmsearch vs the nr database
# the script was run on the cluster
# the output is the rodeo2 output for each segment

for file in LanM_segment*.txt
do
 rodeo2  -out LanM_rodeo_out/rodeo_out_${file%.txt} -j 12 -min 12 -ft 'nucs' -fn 25000 --evaluate_all $file
done
