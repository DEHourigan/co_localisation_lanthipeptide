# made using https://mmseqs.com/latest/userguide.pdf
# PFAM FIRST
# create pfam database for mmseqs this is in pfam dir
    mmseqs convertmsa Pfam-A.full.gz pfam_msa_db 
    mmseqs msa2profile pfam_msa_db pfam_profile --match-mode 1 --threads 10 
    mmseqs createindex pfam_profile tmp -s 7.5 

# make database
    mmseqs createdb ref_seq_proteins_concatenated.faa ref_seq_proteins_concatenated_db
# cluser sequences with mmseqs
    mmseqs cluster --cov-mode 1 --threads 16  ref_seq_proteins_concatenated_db ref_seq_proteins_concatenated.faa tmp
# clusthash redundancy filter sequences with identical length and 100% length overlap
    mmseqs clusthash ref_seq_proteins_concatenated_db  hashDB  --min-seq-id 0.95 --threads 16 # redone on 6th june
    mmseqs clust ref_seq_proteins_concatenated_db hashDB clusterDB # 

# turn cluster results into profile 
    mmseqs createsubdb  clusterDB ref_seq_proteins_concatenated_db clusterDBSeqRepDb        
    mmseqs result2profile clusterDBSeqRepDb ref_seq_proteins_concatenated_db clusterDB refseqclusterProfileDb --threads 24

# get represative of clusterings
    mmseqs profile2repseq refseqclusterProfileDb representative_db --threads 24

# Search now against the created profile database:     loc- /data/san/data0/databases/Pfam/pfam_profile
    mmseqs search representative_db /data/san/data0/databases/Pfam/pfam_profile rep_v_pfam_searchout tmp --threads 12 -k 6 -s 7.5 2>&1 > mmseqs.search.log 
    mmseqs convertalis representative_db /data/san/data0/databases/Pfam/pfam_profile rep_v_pfam_searchout rep_v_pfam_searchout.m8

# to analysis 
mmseqs convertalis ref_seq_proteins_concatenated_db clusterDB hashDb clusterDB.m8
mmseqs convertalis representative_db clusterDB hashDb repDB.m8

#rep_v_pfam_searchout.m8 contains refseq concatenated and PF

# need refseqconcatID and original locus tag

mmseqs createtsv  ref_seq_proteins_concatenated_db representative_db clusterDB  out.tsv
mmseqs createtsv   representative_db clusterDB  out2.tsv  ## pfam can be mapped to whats in this
mmseqs createtsv   representative_db ref_seq_proteins_concatenated_db  out3.tsv