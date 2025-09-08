for i in data/contigs/*.fasta; do
	echo $i
	platon --db /data/san/data0/databases/platon/db --output data/platon_out --threads 4 $i
done
