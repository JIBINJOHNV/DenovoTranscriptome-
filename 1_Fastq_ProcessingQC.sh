##Identify sample names
ls -ltr | awk '{print $9}'  | grep "R1" | sed 's/_R1_001.fastq.gz//g'  > Samplenames.txt 


#Raw data fasq QC
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/Raw_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/${Samplename} ${Samplename}_R1_001.fastq.gz &
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/${Samplename} ${Samplename}_R2_001.fastq.gz &

done < Samplenames.txt 

mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/




##Removing loq quality bases and adaptors using FASTP

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
fastp --in1 ${Samplename}_R1_001.fastq.gz --in2 ${Samplename}_R2_001.fastq.gz \
      --out1 ${Samplename}_Fastp_R1.fastq  --out2 ${Samplename}_Fastp_R2.fastq \
      --length_required 30 --json ${Samplename}_Fastp.json --html ${Samplename}_Fastp.html \
      --correction --detect_adapter_for_pe &
done <Samplenames.txt


mkdir QC/FastpProcessingReports/
mv *html QC/FastpProcessingReports/ ; mv *json QC/FastpProcessingReports/ 

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
gzip ${Samplename}_Fastp_R1.fastq & gzip ${Samplename}_Fastp_R2.fastq 
done <Samplenames.txt



#Qc OF FASTP PROCESSED fasq FILE
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/Fastp_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Fastp_Fastq_QC/${Samplename} ${Samplename}_Fastp_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/Fastp_Fastq_QC/${Samplename} ${Samplename}_Fastp_R2.fastq.gz &

done < Samplenames.txt 

mkdir -p MultiQC/FastQC/Fastp_Fastq_QC
multiqc -o MultiQC/FastQC/Fastp_Fastq_QC/ QC/FastQC/Fastp_Fastq_QC/





##Removing loq quality bases and adaptors using bbduk
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

bbduk.sh in=${Samplename}_Fastp_R1.fastq.gz out=output.fq.gz literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ktrim=r k=20 mink=11 hdist=1
bbduk.sh in=output.fq.gz out=${Samplename}_Fastpbbduk_R1.fastq.gz literal=ATCTCGTATGCCGTCTTCTGCTTG ktrim=r k=15 mink=11 hdist=1 minlength=35

bbduk.sh in=${Samplename}_Fastp_R2.fastq.gz out=output.fq.gz literal=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT ktrim=l k=15 mink=11 hdist=1
bbduk.sh in=output.fq.gz out=${Samplename}_Fastpbbduk_R2.fastq.gz literal=TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTG ktrim=r k=15 mink=11 hdist=1 minlength=35

#trim_galore --length 35  \
#                --paired output_R1.fq.gz output_R2.fq.gz --basename ${Samplename}_TrimgaloreOnlyPair
rm output*
done <Samplenames.txt



#Qc OF FASTP AND BBDUK PROCESSED fasq FILE

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/Fastpbbduk_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Fastpbbduk_Fastq_QC/${Samplename} ${Samplename}_Fastpbbduk_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/Fastpbbduk_Fastq_QC/${Samplename} ${Samplename}_Fastpbbduk_R2.fastq.gz &
done < Samplenames.txt 
wait 

mkdir -p MultiQC/FastQC/Fastpbbduk_Fastq_QC
multiqc -o MultiQC/FastQC/Fastpbbduk_Fastq_QC/ QC/FastQC/Fastpbbduk_Fastq_QC/



###EXTRACTNG qc PASSED PARE END READS

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
fastp --in1 ${Samplename}_Fastpbbduk_R1.fastq.gz --in2 ${Samplename}_Fastpbbduk_R2.fastq.gz \
      --out1 ${Samplename}_Paired_Fastpbbduk_R1.fastq.gz  --out2 ${Samplename}_Paired_Fastpbbduk_R2.fastq.gz \
      --length_required 30 --json ${Samplename}_Paired_Fastpbbduk.json --html ${Samplename}_Paired_Fastpbbduk.html \
      --correction --disable_adapter_trimming &
done < Samplenames.txt 



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/Fastpbbduk_PAIRED_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Fastpbbduk_PAIRED_Fastq_QC/${Samplename} ${Samplename}_Paired_Fastpbbduk_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/Fastpbbduk_PAIRED_Fastq_QC/${Samplename} ${Samplename}_Paired_Fastpbbduk_R2.fastq.gz &
done < Samplenames.txt 
wait 

mkdir -p MultiQC/FastQC/Fastpbbduk_PAIRED_Fastq_QC
multiqc -o MultiQC/FastQC/Fastpbbduk_PAIRED_Fastq_QC/ QC/FastQC/Fastpbbduk_PAIRED_Fastq_QC/


mkdir QC/Fastpbbduk_PAIREDProcessingReports/
mv *html QC/Fastpbbduk_PAIREDProcessingReports/ ; mv *json QC/Fastpbbduk_PAIREDProcessingReports/ 



mkdir -p ProcessedFastQ_Files/Paired_Fastpbbduk
mv *_Paired_Fastpbbduk_* ProcessedFastQ_Files/Paired_Fastpbbduk

mkdir -p ProcessedFastQ_Files/Fastpbbduk
mv *_Fastpbbduk_* ProcessedFastQ_Files/Fastpbbduk

mkdir -p ProcessedFastQ_Files/Fastp
mv *_Fastp_* ProcessedFastQ_Files/Fastp


####################################-------------------------------------------------------------###################################################################



##MERING FASTQ FILES
zcat C1_S1_L003_Paired_Fastpbbduk_R1.fastq.gz C2_S2_L003_Paired_Fastpbbduk_R1.fastq.gz C2_S11_L002_Paired_Fastpbbduk_R1.fastq.gz | gzip > Control-IP-TPaired_Fastpbbduk_R1.fastq.gz & 
zcat C1_S1_L003_Paired_Fastpbbduk_R2.fastq.gz C2_S2_L003_Paired_Fastpbbduk_R2.fastq.gz C2_S11_L002_Paired_Fastpbbduk_R2.fastq.gz | gzip > Control-IP-TPaired_Fastpbbduk_R2.fastq.gz &

zcat T1_S3_L003_Paired_Fastpbbduk_R1.fastq.gz T2_S4_L003_Paired_Fastpbbduk_R1.fastq.gz | gzip > Treated-IP-TPaired_Fastpbbduk_R1.fastq.gz & 
zcat T1_S3_L003_Paired_Fastpbbduk_R2.fastq.gz T2_S4_L003_Paired_Fastpbbduk_R2.fastq.gz | gzip > Treated-IP-TPaired_Fastpbbduk_R2.fastq.gz & 

zcat input-C1_S6_L003_Paired_Fastpbbduk_R1.fastq.gz input-C2_S7_L003_Paired_Fastpbbduk_R1.fastq.gz input-C1_S12_L002_Paired_Fastpbbduk_R1.fastq.gz input-C2_S13_L002_Paired_Fastpbbduk_R1.fastq.gz | gzip > Input-Control-TPaired_Fastpbbduk_R1.fastq.gz &
zcat input-C1_S6_L003_Paired_Fastpbbduk_R2.fastq.gz input-C2_S7_L003_Paired_Fastpbbduk_R2.fastq.gz input-C1_S12_L002_Paired_Fastpbbduk_R2.fastq.gz input-C2_S13_L002_Paired_Fastpbbduk_R2.fastq.gz | gzip > Input-Control-TPaired_Fastpbbduk_R2.fastq.gz &

zcat input-T1_S8_L003_Paired_Fastpbbduk_R1.fastq.gz  input-T2_S9_L003_Paired_Fastpbbduk_R1.fastq.gz input-T1_S14_L002_Paired_Fastpbbduk_R1.fastq.gz input-T2_S15_L002_Paired_Fastpbbduk_R1.fastq.gz | gzip > Input-Treated-TPaired_Fastpbbduk_R1.fastq.gz &
zcat input-T1_S8_L003_Paired_Fastpbbduk_R2.fastq.gz  input-T2_S9_L003_Paired_Fastpbbduk_R2.fastq.gz input-T1_S14_L002_Paired_Fastpbbduk_R2.fastq.gz input-T2_S15_L002_Paired_Fastpbbduk_R2.fastq.gz | gzip > Input-Treated-TPaired_Fastpbbduk_R2.fastq.gz &



mkdir  -p MergedFastQFiles/Paired_Fastpbbduk
mv *-TPaired_Fastpbbduk* MergedFastQFiles/Paired_Fastpbbduk


ls -ltr MergedFastQFiles/Paired_Fastpbbduk/ | awk '{print $9}' |grep "R1"| sed 's/-TPaired_Fastpbbduk_R1.fastq.gz//g' > MergedSampleNames.txt



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/merged_Fastpbbduk_PAIRED/${Samplename}
fastqc -t 5 -o QC/FastQC/merged_Fastpbbduk_PAIRED/${Samplename} MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/merged_Fastpbbduk_PAIRED/${Samplename} MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R2.fastq.gz &
done < MergedSampleNames.txt 
wait 


mkdir -p MultiQC/FastQC/merged_Fastpbbduk_PAIRED
multiqc -o MultiQC/FastQC/merged_Fastpbbduk_PAIRED/ QC/FastQC/merged_Fastpbbduk_PAIRED/







#####################################---------------------------------------------------------------------------###################################################

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
fastp --in1 MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R1.fastq.gz --in2 MergedFastQFiles/Paired_Fastpbbduk/${Samplename}-TPaired_Fastpbbduk_R2.fastq.gz \
      --out1 ${Samplename}_QCPaired_Fastpbbduk_R1.fastq.gz  --out2 ${Samplename}_QCPaired_Fastpbbduk_R2.fastq.gz \
      --length_required 30 --json ${Samplename}_QCPaired_Fastpbbduk.json --html ${Samplename}_QCPaired_Fastpbbduk.html \
      --correction --disable_adapter_trimming --cut_front --cut_tail --cut_window_size 3 --cut_mean_quality 15 &
done < MergedSampleNames.txt



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/merged_QCPaired_Fastpbbduk/${Samplename}
fastqc -t 5 -o QC/FastQC/merged_QCPaired_Fastpbbduk/${Samplename} ${Samplename}_QCPaired_Fastpbbduk_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/merged_QCPaired_Fastpbbduk/${Samplename} ${Samplename}_QCPaired_Fastpbbduk_R2.fastq.gz &
done < MergedSampleNames.txt 






























































########################################---------------------------------------------------------------------------##################################################






while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
gzip ${Samplename}_Fastp_R1.fastq & gzip ${Samplename}_Fastp_R2.fastq 
done <Samplenames.txt




#Fastprocessed fastq file QC
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/FastP_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/FastP_Fastq_QC/${Samplename} ${Samplename}_Fastp_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/FastP_Fastq_QC/${Samplename} ${Samplename}_Fastp_R2.fastq.gz &
done < Samplenames.txt 
wait 




mkdir -p MultiQC/FastQC/FastP_Fastq_QC
multiqc -o MultiQC/FastQC/FastP_Fastq_QC/ QC/FastQC/FastP_Fastq_QC/








#######################Merge


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/Trimgalore/${Samplename}
fastqc -t 5 -o QC/FastQC/Trimgalore/${Samplename} ${Samplename}_Paired_Fastpbbduk_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/Trimgalore/${Samplename} ${Samplename}_rimgaloreOnlyPair_val_2.fq.gz &
done < SampleName.txt 


mkdir -p MultiQC/FastQC/Trimgalore_QC
multiqc -o MultiQC/FastQC/Trimgalore_QC/ QC/FastQC/Trimgalore/


mkdir -p ProcessedFastq/Trimgalore
mv *rimgaloreOnlyPair_val* ProcessedFastq/Trimgalore



##MERING FASTQ FILES
zcat C1_S1_L003_Paired_Fastpbbduk_R1.fastq.gz C2_S2_L003_rimgaloreOnlyPair_val_1.fq.gz C2_S11_L002_rimgaloreOnlyPair_val_1.fq.gz | gzip > Control-IP-TrimgaloreOnlyPair_val_1.fq.gz & 
zcat C1_S1_L003_rimgaloreOnlyPair_val_2.fq.gz C2_S2_L003_rimgaloreOnlyPair_val_2.fq.gz C2_S11_L002_rimgaloreOnlyPair_val_2.fq.gz | gzip > Control-IP-TrimgaloreOnlyPair_val_2.fq.gz &

zcat T1_S3_L003_rimgaloreOnlyPair_val_1.fq.gz T2_S4_L003_rimgaloreOnlyPair_val_1.fq.gz | gzip > Treated-IP-TrimgaloreOnlyPair_val_1.fq.gz & 
zcat T1_S3_L003_rimgaloreOnlyPair_val_2.fq.gz T2_S4_L003_rimgaloreOnlyPair_val_2.fq.gz | gzip > Treated-IP-TrimgaloreOnlyPair_val_2.fq.gz & 

zcat input-C1_S6_L003_rimgaloreOnlyPair_val_1.fq.gz input-C2_S7_L003_rimgaloreOnlyPair_val_1.fq.gz input-C1_S12_L002_rimgaloreOnlyPair_val_1.fq.gz input-C2_S13_L002_rimgaloreOnlyPair_val_1.fq.gz | gzip > Input-Control-TrimgaloreOnlyPair_val_1.fq.gz &
zcat input-C1_S6_L003_rimgaloreOnlyPair_val_2.fq.gz input-C2_S7_L003_rimgaloreOnlyPair_val_2.fq.gz input-C1_S12_L002_rimgaloreOnlyPair_val_2.fq.gz input-C2_S13_L002_rimgaloreOnlyPair_val_2.fq.gz | gzip > Input-Control-TrimgaloreOnlyPair_val_2.fq.gz &

zcat input-T1_S8_L003_rimgaloreOnlyPair_val_1.fq.gz  input-T2_S9_L003_rimgaloreOnlyPair_val_1.fq.gz input-T1_S14_L002_rimgaloreOnlyPair_val_1.fq.gz input-T2_S15_L002_rimgaloreOnlyPair_val_1.fq.gz | gzip > Input-Treated-TrimgaloreOnlyPair_val_1.fq.gz &
zcat input-T1_S8_L003_rimgaloreOnlyPair_val_2.fq.gz  input-T2_S9_L003_rimgaloreOnlyPair_val_2.fq.gz input-T1_S14_L002_rimgaloreOnlyPair_val_2.fq.gz input-T2_S15_L002_rimgaloreOnlyPair_val_2.fq.gz | gzip > Input-Treated-TrimgaloreOnlyPair_val_2.fq.gz &


mkdir -p ~/Downloads/ProcessedFastq/MergedFastq
mv *-TrimgaloreOnlyPair_val_*.fq.gz ~/Downloads/ProcessedFastq/MergedFastq


ls  -ltr *-TrimgaloreOnlyPair_val_1.fq.gz | awk '{print $9}' | sed 's/-TrimgaloreOnlyPair_val_1.fq.gz//g' > SampleNames.txt

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
mkdir -p QC/FastQC/Merged_Fastq/${Samplename}
fastqc -t 5 -o QC/FastQC/Merged_Fastq/${Samplename} ${Samplename}-TrimgaloreOnlyPair_val_1.fq.gz &
fastqc -t 5 -o QC/FastQC/Merged_Fastq/${Samplename} ${Samplename}-TrimgaloreOnlyPair_val_2.fq.gz &
done < SampleNames.txt



























ls -ltr *001.fastq.gz   | sort |awk '{print $9}'  | sort | sed 's/_R1_001.fastq.gz//g' | sed 's/_R2_001.fastq.gz//g' | uniq > SampleName.txt




while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

echo "Step1: Prealignment QC" 

mkdir -p QC/FastQC/Raw_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/${Samplename} ${Samplename}_R1_001.fastq.gz &
fastqc -t 5 -o QC/FastQC/Raw_Fastq_QC/${Samplename} ${Samplename}_R2_001.fastq.gz &

done < SampleName.txt 
wait 


mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/




while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
java -jar  Trimmomatic-0.39/trimmomatic-0.39.jar PE ${Samplename}_R1_001.fastq.gz ${Samplename}_R2_001.fastq.gz  \
           ${Samplename}_Paired_R1.fastq.gz ${Samplename}_UnPaired_R1.fastq.gz \
           ${Samplename}_Paired_R2.fastq.gz ${Samplename}_UnPaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:3:15 LEADING:15 TRAILING:15 MINLEN:30  
done < SampleName.txt 
wait 

java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE input-C2_S13_L002_R1_001.fastq.gz input-C2_S13_L002_R1_001.fastq.gz \
           input-C2_Paired_R2.fastq.gz input-C2_UnPaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:3:15 LEADING:15 TRAILING:15 MINLEN:30  


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
trim_galore --length 35  \
                --paired ${Samplename}_R1_001.fastq.gz ${Samplename}_R2_001.fastq.gz --basename ${Samplename}_trimgaloreOnlyPair &
done < SampleName.txt 
wait 



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
fastp --in1  ${Samplename}_trimgaloreOnlyPair_val_1.fq.gz   --in2  ${Samplename}_trimgaloreOnlyPair_val_2.fq.gz  \
      --out1 ${Samplename}_QC_R1.fastq.gz           --out2 ${Samplename}_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --json ${Samplename}_QC.json --html ${Samplename}_QC.html -A --correction --length_required 30  &
done < SampleName.txt 
wait 


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/FastP_Fastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/FastP_Fastq_QC/${Samplename} ${Samplename}_QC_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/FastP_Fastq_QC/${Samplename} ${Samplename}_QC_R2.fastq.gz &

done < SampleName.txt 





mkdir -p ProcessedFASTQ/trimgalore
mv *trimgaloreOnlyPair*fq.gz ProcessedFASTQ/trimgalore

mkdir -p QC/FastqProcessing/Trimgalore
mv *fastq.gz_trimming_report.txt QC/FastqProcessing/Trimgalore




##MERING FASTQ FILES
zcat C1_S1_L003_QC_R1.fastq.gz C2_S2_L003_QC_R1.fastq.gz C2_S11_L002_QC_R1.fastq.gz | gzip > Control-IP-QC_R1.fastq.gz & 
zcat C1_S1_L003_QC_R2.fastq.gz C2_S2_L003_QC_R2.fastq.gz C2_S11_L002_QC_R2.fastq.gz | gzip > Control-IP-QC_R2.fastq.gz &

zcat T1_S3_L003_QC_R1.fastq.gz T2_S4_L003_QC_R1.fastq.gz | gzip > Treated-IP-QC_R1.fastq.gz & 
zcat T1_S3_L003_QC_R2.fastq.gz T2_S4_L003_QC_R2.fastq.gz | gzip > Treated-IP-QC_R2.fastq.gz & 

zcat input-C1_S6_L003_QC_R1.fastq.gz input-C2_S7_L003_QC_R1.fastq.gz input-C1_S12_L002_QC_R1.fastq.gz input-C2_S13_L002_QC_R1.fastq.gz | gzip > Input-Control-QC_R1.fastq.gz &
zcat input-C1_S6_L003_QC_R2.fastq.gz input-C2_S7_L003_QC_R2.fastq.gz input-C1_S12_L002_QC_R2.fastq.gz input-C2_S13_L002_QC_R2.fastq.gz | gzip > Input-Control-QC_R2.fastq.gz &

zcat input-T1_S8_L003_QC_R1.fastq.gz  input-T2_S9_L003_QC_R1.fastq.gz input-T1_S14_L002_QC_R1.fastq.gz input-T2_S15_L002_QC_R1.fastq.gz | gzip > Input-Treated-QC_R1.fastq.gz &
zcat input-T1_S8_L003_QC_R2.fastq.gz  input-T2_S9_L003_QC_R2.fastq.gz input-T1_S14_L002_QC_R2.fastq.gz input-T2_S15_L002_QC_R2.fastq.gz | gzip > Input-Treated-QC_R2.fastq.gz &




ls -ltr *QC_R1.fastq.gz  | awk '{print $9}' | sed 's/_R1.fastq.gz//g' > SampleName2.txt


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/FinalMergedFastq_QC/${Samplename}
fastqc -t 5 -o QC/FastQC/FinalMergedFastq_QC/${Samplename} ${Samplename}_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/FinalMergedFastq_QC/${Samplename} ${Samplename}_R2.fastq.gz &

done < SampleName2.txt 


mkdir -p ProcessedFASTQ/Fastp
mv *_L00*_QC_R*.fastq.gz ProcessedFASTQ/Fastp

mkdir -p QC/FastqProcessing/Fastp
mv *.html QC/FastqProcessing/Fastp ; mv *.json QC/FastqProcessing/Fastp



##MERING FASTQ FILES





























































zcat input-C1_S6_L003_trimgalore_R1.fastq.gz input-C2_S7_L003_trimgalore_R1.fastq.gz input-C1_S12_L002_trimgalore_R1.fastq.gz input-C2_S13_L002_trimgalore_R1.fastq.gz | gzip > Input-Control_trimgalore_R1.fastq.gz & 
zcat input-C1_S6_L003_trimgalore_R2.fastq.gz input-C2_S7_L003_trimgalore_R2.fastq.gz input-C1_S12_L002_trimgalore_R2.fastq.gz  input-C2_S13_L002_trimgalore_R2.fastq.gz | gzip > Input-Control_trimgalore_R2.fastq.gz & 

zcat input-T1_S14_L002_trimgalore_R1.fastq.gz input-T1_S8_L003_trimgalore_R1.fastq.gz input-T2_S15_L002_trimgalore_R1.fastq.gz input-T2_S9_L003_trimgalore_R1.fastq.gz | gzip >Input-Treated_trimgalore_R1.fastq.gz & 
zcat input-T1_S14_L002_trimgalore_R2.fastq.gz input-T1_S8_L003_trimgalore_R2.fastq.gz input-T2_S15_L002_trimgalore_R2.fastq.gz input-T2_S9_L003_trimgalore_R2.fastq.gz | gzip >Input-Treated_trimgalore_R2.fastq.gz &

zcat C1_S1_L003_trimgalore_R1.fastq.gz C2_S11_L002_trimgalore_R1.fastq.gz C2_S2_L003_trimgalore_R1.fastq.gz | gzip > Control_IP_trimgalore_R1.fastq.gz &
zcat C1_S1_L003_trimgalore_R2.fastq.gz C2_S11_L002_trimgalore_R2.fastq.gz C2_S2_L003_trimgalore_R2.fastq.gz | gzip > Control_IP_trimgalore_R2.fastq.gz &

zcat T1_S3_L003_trimgalore_R1.fastq.gz T2_S4_L003_trimgalore_R1.fastq.gz | gzip >Treated_IP_trimgalore_R1.fastq.gz &
zcat T1_S3_L003_trimgalore_R2.fastq.gz T2_S4_L003_trimgalore_R2.fastq.gz | gzip >Treated_IP_trimgalore_R2.fastq.gz &


mkdir -p QC/FastQProcessing/Trimgalore
mv *trimming_* QC/FastQProcessing/Trimgalore
mv *_QC.json QC/FastQProcessing/Trimgalore
mv *_QC.html QC/FastQProcessing/Trimgalore


ls -ltr *trimgalore_R1.fastq.gz  | sed 's/_trimgalore_R1.fastq.gz//g' | awk '{print $9}' > SampleName.txt






while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
trim_galore --length 35 --retain_unpaired --length_1 50 --length_2 50 \
                --paired ${Samplename}_R1_001.fastq.gz ${Samplename}_R2_001.fastq.gz --basename ${Samplename}_trimgalore &
done < SampleName.txt 
wait 

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
zcat ${Samplename}_trimgalore_val_1.fq.gz ${Samplename}_trimgalore_R1_unpaired_1.fq.gz | gzip >${Samplename}_trimgalore_R1.fastq.gz &
zcat ${Samplename}_trimgalore_val_2.fq.gz ${Samplename}_trimgalore_R2_unpaired_2.fq.gz | gzip >${Samplename}_trimgalore_R2.fastq.gz & 

done < SampleName.txt 
wait 



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 
fastp --in1  ${Samplename}_trimgalore_R1.fastq.gz   --in2  ${Samplename}_trimgalore_R2.fastq.gz  \
      --out1 ${Samplename}_QC_R1.fastq.gz           --out2 ${Samplename}_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --json ${Samplename}_QC.json --html ${Samplename}_QC.html -A --correction --length_required 30 
done < SampleName.txt 
wait 



mkdir -p ProcessedFASTQ/fastp
<pre>mv *_QC_R[12].fastq.gz ProcessedFASTQ/fastp</pre>




##MERING FASTQ FILES
zcat C1_S1_L003_QC_R1.fastq.gz C2_S2_L003_QC_R1.fastq.gz C2_S11_L002_QC_R1.fastq.gz | gzip > Control_IP_QC_R1.fastq.gz & 
zcat C1_S1_L003_QC_R2.fastq.gz C2_S2_L003_QC_R2.fastq.gz C2_S11_L002_QC_R2.fastq.gz | gzip > Control_IP_QC_R2.fastq.gz &

zcat T1_S3_L003_QC_R1.fastq.gz T2_S4_L003_QC_R1.fastq.gz | gzip > Treated_IP_QC_R1.fastq.gz & 
zcat T1_S3_L003_QC_R2.fastq.gz T2_S4_L003_QC_R2.fastq.gz | gzip > Treated_IP_QC_R2.fastq.gz & 

zcat input-C1_S6_L003_QC_R1.fastq.gz input-C2_S7_L003_QC_R1.fastq.gz input-C1_S12_L002_QC_R1.fastq.gz input-C2_S13_L002_QC_R1.fastq.gz | gzip > Input_Control_QC_R1.fastq.gz &
zcat input-C1_S6_L003_QC_R2.fastq.gz input-C2_S7_L003_QC_R2.fastq.gz input-C1_S12_L002_QC_R2.fastq.gz input-C2_S13_L002_QC_R2.fastq.gz | gzip > Input_Control_QC_R2.fastq.gz &

zcat input-T1_S8_L003_QC_R1.fastq.gz  input-T2_S9_L003_QC_R1.fastq.gz input-T1_S14_L002_QC_R1.fastq.gz input-T2_S15_L002_QC_R1.fastq.gz | gzip > Input_Treated_QC_R1.fastq.gz &
zcat input-T1_S8_L003_QC_R2.fastq.gz  input-T2_S9_L003_QC_R2.fastq.gz input-T1_S14_L002_QC_R2.fastq.gz input-T2_S15_L002_QC_R2.fastq.gz | gzip > Input_Treated_QC_R2.fastq.gz &




ls *gz | egrep -v "001|002" | grep "R1.fastq.gz" | sed 's/_R1.fastq.gz//g'  > Sample2.txt

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait 

mkdir -p QC/FastQC/QC_Passed/${Samplename}
fastqc -t 5 -o QC/FastQC/QC_Passed/${Samplename} ${Samplename}_R1.fastq.gz &
fastqc -t 5 -o QC/FastQC/QC_Passed/${Samplename} ${Samplename}_R2.fastq.gz &

done < Sample2.txt
wait 


#MultiQC
mkdir -p QC/MultiQC/QC_Passed/

multiqc --outdir QC/MultiQC/QC_Passed/  QC/FastQC/QC_Passed/ 







##MERING FASTQ FILES
zcat C1_S1_L003_R1_001.fastq.gz C2_S2_L003_R1_001.fastq.gz C2_S11_L002_R1_001.fastq.gz | gzip > Control_IP_R1.fastq.gz
zcat C1_S1_L003_R2_001.fastq.gz C2_S2_L003_R2_001.fastq.gz C2_S11_L002_R2_001.fastq.gz | gzip > Control_IP_R2.fastq.gz

zcat T1_S3_L003_R1_001.fastq.gz T2_S4_L003_R1_001.fastq.gz | gzip > Treated_IP_R1.fastq.gz
zcat T1_S3_L003_R2_001.fastq.gz T2_S4_L003_R2_001.fastq.gz | gzip > Treated_IP_R2.fastq.gz

zcat input-C1_S6_L003_R1_001.fastq.gz input-C2_S7_L003_R1_001.fastq.gz input-C1_S12_L002_R1_001.fastq.gz input-C2_S13_L002_R1_001.fastq.gz | gzip > Input_Control_R1.fastq.gz
zcat input-C1_S6_L003_R2_001.fastq.gz input-C2_S7_L003_R2_001.fastq.gz input-C1_S12_L002_R2_001.fastq.gz input-C2_S13_L002_R2_001.fastq.gz | gzip > Input_Control_R2.fastq.gz

zcat input-T1_S8_L003_R1_001.fastq.gz  input-T2_S9_L003_R1_001.fastq.gz input-T1_S14_L002_R1_001.fastq.gz input-T2_S15_L002_R1_001.fastq.gz | gzip > Input_Treated_R1.fastq.gz
zcat input-T1_S8_L003_R2_001.fastq.gz  input-T2_S9_L003_R2_001.fastq.gz input-T1_S14_L002_R2_001.fastq.gz input-T2_S15_L002_R2_001.fastq.gz | gzip > Input_Treated_R2.fastq.gz



#Raw fastq QC
mkdir -p QC/FastQC/RawQC/

fastqc -o QC/FastQC/RawQC/ Treated_IP_R1.fastq.gz Treated_IP_R2.fastq.gz
fastqc -o QC/FastQC/RawQC/ Control_IP_R1.fastq.gz Control_IP_R2.fastq.gz

fastqc -o QC/FastQC/RawQC/ Input_Control_R1.fastq.gz Input_Control_R2.fastq.gz
fastqc -o QC/FastQC/RawQC/ Input_Treated_R1.fastq.gz Input_Treated_R2.fastq.gz



##Trimming using trim_galore
trim_galore --length 35 --retain_unpaired --length_1 50 --length_2 50 --paired Treated_IP_R1.fastq.gz     Treated_IP_R2.fastq.gz    --fastqc --basename Treated_IP_trimgalore 
trim_galore --length 35 --retain_unpaired --length_1 50 --length_2 50 --paired Control_IP_R1.fastq.gz     Control_IP_R2.fastq.gz    --fastqc --basename Control_IP_trimgalore 
trim_galore --length 35 --retain_unpaired --length_1 50 --length_2 50 --paired Input_Control_R1.fastq.gz  Input_Control_R2.fastq.gz --fastqc --basename Input_Control_trimgalore 
trim_galore --length 35 --retain_unpaired --length_1 50 --length_2 50 --paired Input_Treated_R1.fastq.gz  Input_Treated_R2.fastq.gz --fastqc --basename Input_Treated_trimgalore 




/home/jibin/Softwar/atria-3.2.1/bin/atria --detect-adapter \
--read1 Treated_IP_R1.fastq.gz --read2 Treated_IP_R2.fastq.gz --output-dir Treated_IP --quality-kmer 3 --length-range 30:500




zcat Input_Control_trimgalore_val_1.fq.gz Input_Control_trimgalore_R1_unpaired_1.fq.gz | gzip >Input_Control_trimgalore_R1.fastq.gz 
zcat Input_Control_trimgalore_val_2.fq.gz Input_Control_trimgalore_R2_unpaired_2.fq.gz | gzip >Input_Control_trimgalore_R2.fastq.gz 

zcat Control_IP_trimgalore_val_1.fq.gz Control_IP_trimgalore_R1_unpaired_1.fq.gz | gzip > Control_IP_trimgalore_R1.fastq.gz 
zcat Control_IP_trimgalore_val_2.fq.gz Control_IP_trimgalore_R2_unpaired_2.fq.gz | gzip > Control_IP_trimgalore_R2.fastq.gz 

zcat Input_Treated_trimgalore_val_1.fq.gz Input_Treated_trimgalore_R1_unpaired_1.fq.gz | gzip > Input_Treated_trimgalore_R1.fastq.gz 
zcat Input_Treated_trimgalore_val_2.fq.gz Input_Treated_trimgalore_R2_unpaired_2.fq.gz | gzip > Input_Treated_trimgalore_R2.fastq.gz 

zcat Treated_IP_trimgalore_val_1.fq.gz Treated_IP_trimgalore_R1_unpaired_1.fq.gz | gzip > Treated_IP_trimgalore_R1.fastq.gz 
zcat Treated_IP_trimgalore_val_2.fq.gz Treated_IP_trimgalore_R2_unpaired_2.fq.gz | gzip > Treated_IP_trimgalore_R2.fastq.gz 


##Trimming using trim_galore
fastp --in1  Input_Control_trimgalore_R1.fastq.gz   --in2  Input_Control_trimgalore_R2.fastq.gz  \
      --out1 Input_Control_QC_R1.fastq.gz           --out2 Input_Control_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --json Input_Control_QC.json --html Input_Control_QC.html -A --correction --length_required 30 

fastp --in1  Control_IP_trimgalore_R1.fastq.gz    --in2  Control_IP_trimgalore_R2.fastq.gz \
      --out1 Control_IP_QC_R1.fastq.gz --out2 Control_IP_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --html Control_IP_QC.html --json Control_IP_QC.json -A --correction --length_required 30  

fastp --in1  Input_Treated_trimgalore_R1.fastq.gz    --in2  Input_Treated_trimgalore_R2.fastq.gz \
      --out1 Input_Treated_QC_R1.fastq.gz --out2 Input_Treated_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --html Input_Treated_QC.html --json Input_Treated_QC.json -A --correction --length_required 30  

fastp --in1   Treated_IP_trimgalore_R1.fastq.gz    --in2   Treated_IP_trimgalore_R2.fastq.gz \
      --out1  Treated_IP_QC_R1.fastq.gz --out2  Treated_IP_QC_R2.fastq.gz \
      --cut_front --cut_tail --cut_window_size 3 --html  Treated_IP_QC.html --json  Treated_IP_QC.json -A --correction --length_required 30 



mkdir -p ProcessedFastQ/trimgalore
mv *trimgalore*gz  ProcessedFastQ/trimgalore/

mkdir -p QC/fastp ; mv *json  QC/fastp/ ; mv *html  QC/fastp/



mkdir -p QC/FastQC/FastptrimgaloreQC/

fastqc -o QC/FastQC/FastptrimgaloreQC/ Treated_IP_QC_R1.fastq.gz    Treated_IP_QC_R2.fastq.gz
fastqc -o QC/FastQC/FastptrimgaloreQC/ Control_IP_QC_R1.fastq.gz    Control_IP_QC_R2.fastq.gz
fastqc -o QC/FastQC/FastptrimgaloreQC/ Input_Control_QC_R1.fastq.gz Input_Control_QC_R2.fastq.gz
fastqc -o QC/FastQC/FastptrimgaloreQC/ Input_Treated_QC_R1.fastq.gz Input_Treated_QC_R2.fastq.gz


mkdir -p MultiQC/FastQC/FastptrimgaloreQC/  
multiqc --outdir MultiQC/FastQC/FastptrimgaloreQC/  QC/FastQC/FastptrimgaloreQC/ 

################
