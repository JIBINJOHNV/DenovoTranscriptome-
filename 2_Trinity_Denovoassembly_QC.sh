


docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
      --seqType fq \
      --left `pwd`/MergedFastQFiles/Paired_Fastpbbduk/Control-IP-TPaired_Fastpbbduk_R1.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Input-Control-TPaired_Fastpbbduk_R1.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Input-Treated-TPaired_Fastpbbduk_R1.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Treated-IP-TPaired_Fastpbbduk_R1.fastq.gz \
      --right `pwd`/MergedFastQFiles/Paired_Fastpbbduk/Control-IP-TPaired_Fastpbbduk_R2.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Input-Control-TPaired_Fastpbbduk_R2.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Input-Treated-TPaired_Fastpbbduk_R2.fastq.gz,`pwd`/MergedFastQFiles/Paired_Fastpbbduk/Treated-IP-TPaired_Fastpbbduk_R2.fastq.gz \
      --max_memory 5G --CPU 4 --output `pwd`/Trinity &


##QUANTFY 
mkdir -p Quantification/kallisto

docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq \
 /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts `pwd`/Trinity.Trinity.fasta \
--est_method kallisto \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference


Input-Control-Tpaired_Fastpbbduk_R1.fastq.gz
Input-Control-TPaired_Fastpbbduk_R1.fastq.gz


mkdir -p  Quantification/kallisto/Treated_1/
mkdir -p  Quantification/kallisto/Control_1/


docker run -it -v`pwd`:`pwd` -v`pwd`:/Data/ trinityrnaseq/trinityrnaseq 


/usr/local/bin/util/align_and_estimate_abundance.pl \
    --transcripts /Data/Trinity.Trinity.fasta \
    --samples_file /Data/SampleDetails.txt --seqTyp fq \
     --est_method kallisto --aln_method bowtie2 --trinity_mode --output_dir /Data/Quantification/kallisto/


/usr/local/bin/util/support_scripts/kallisto_trans_to_gene_results.pl \
Treated_1/abundance.tsv /Data/Trinity.Trinity.fasta.gene_trans_map > /Data/Quantification/kallisto/Treated_1/abundance.tsv.genes

/usr/local/bin/util/support_scripts/kallisto_trans_to_gene_results.pl \
Control_1/abundance.tsv /Data/Trinity.Trinity.fasta.gene_trans_map > /Data/Quantification/kallisto/Control_1/abundance.tsv.genes





############


docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq \
    /usr/local/bin/bowtie2-build `pwd`/Trinity.Trinity.fasta \
     `pwd`/Trinity.Trinity.fasta

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/run_bowtie2.pl \
                --target `pwd`/Trinity.Trinity.fasta  \
                --left `pwd`/MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R1.fastq.gz --right `pwd`/MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R2.fastq.gz  \
                | samtools view -Sb - | samtools sort - -o `pwd`/${Samplename}_coordSorted.bam & 
done < MergedSampleNames.txt



docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.3_cv1  busco -m transcriptome -i /busco_wd/Trinity.Trinity.fasta \
  -o /busco_wd/Trinity.Trinity.fasta.buscoviridi  -l viridiplantae_odb10  -f


##Counting Full Length Trinity Transcripts
# Build a blastable database:


docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq makeblastdb -in `pwd`/Database/uniprot_sprot.fasta -dbtype prot

# Perform the blast search, reporting only the top alignment:
docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq blastx \
            -query trinity_out_dir.fasta \
            -db `pwd`/Database/uniprot_sprot.fasta -out `pwd`/blastx.outfmt6 \
            -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6

$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 Trinity.fasta uniprot_sprot.fasta



#Grouping Blast Hits to Improve Sequence Coverage
