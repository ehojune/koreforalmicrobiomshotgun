#!/bin/bash


### paths ###

rawdata="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/raw_data/"
scripts="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/scripts/"
output="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data"
hg38="GRCh38_noalt_as"
#samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2"
samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2 KorefMicrobiomeShotgun3"

### programs ###

fastp="/BiO/Access/ehojune/anaconda3/bin/fastp"
megahit="/BiO/Access/ehojune/anaconda3/bin/megahit"
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py"
quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"
bowtie2="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bowtie2-2.5.1/bowtie2"
samtools="/BiO/Access/ehojune/anaconda3/bin/samtools"
bedtools="/BiO/Access/ehojune/anaconda3/bin/bedtools"


### commands ###


:<<'END'

#1. fastp before setting trimming criteria

output_0=${output}/0_fastp_beforeTrimming/
mkdir -p $output_0

echo "start fastp before"
for sample in $samples
do
  mkdir -p ${output_0}${sample}
  $fastp -i ${rawdata}${sample}/${sample}_1.fq.gz -I ${rawdata}${sample}/${sample}_2.fq.gz -o ${output_0}${sample}/${sample}_1.fq.gz -O ${output_0}${sample}/${sample}_2.fq.gz --html ${output_0}${sample}/${sample}.html --json ${output_0}${sample}/${sample}.json --report_title ${sample}_before  --overrepresentation_analysis --thread 16
done
echo "done"




#2. fastp after setting trimming criteria

output_1=${output}/1_fastp_afterTrimming/
mkdir -p $output_1

for sample in $samples
do
  mkdir -p ${output_1}${sample}
done




echo "start fastp after"

touch ${scripts}KorefMicrobiomeShotgun1_fastp_after.sh
echo $fastp -i ${rawdata}${sample1}/${sample1}_1.fq.gz -I ${rawdata}${sample1}/${sample1}_2.fq.gz -o ${output_1}${sample1}/${sample1}_1.fq.gz -O ${output_1}${sample1}/${sample1}_2.fq.gz --html ${output_1}${sample1}/${sample1}.html --json ${output_1}${sample1}/${sample1}.json --report_title ${sample1}_after  --overrepresentation_analysis --trim_front1 18 --trim_tail1 5 --trim_front2 18 --trim_tail2 3 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20 > ${scripts}KorefMicrobiomeShotgun1_fastp_after.sh 
qsub -pe smp 16 -cwd ${scripts}KorefMicrobiomeShotgun1_fastp_after.sh








touch ${scripts}KorefMicrobiomeShotgun2_fastp_after.sh
echo $fastp -i ${rawdata}${sample2}/${sample2}_1.fq.gz -I ${rawdata}${sample2}/${sample2}_2.fq.gz -o ${output_1}${sample2}/${sample2}_1.fq.gz -O ${output_1}${sample2}/${sample2}_2.fq.gz --html ${output_1}${sample2}/${sample2}.html --json ${output_1}${sample2}/${sample2}.json --report_title ${sample2}_after  --overrepresentation_analysis --trim_front1 10 --trim_front2 10 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20 > ${scripts}KorefMicrobiomeShotgun2_fastp_after.sh
qsub -pe smp 16 -cwd ${scripts}KorefMicrobiomeShotgun2_fastp_after.sh









mkdir -p ${output_1}KorefMicrobiomeShotgun3
cat ${output_1}${sample1}/${sample1}_1.fq.gz ${output_1}${sample2}/${sample2}_1.fq.gz > ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_1.fq.gz
cat ${output_1}${sample1}/${sample1}_2.fq.gz ${output_1}${sample2}/${sample2}_2.fq.gz > ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_2.fq.gz





output_1=${output}/1_fastp_afterTrimming/

touch ${scripts}KorefMicrobiomeShotgun3_fastp.sh
echo $fastp -i ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_1.fq.gz -I ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_2.fq.gz -o ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_1.fq.gz -O ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_2.fq.gz --html ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.html --json ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.json --report_title KorefMicrobiomeShotgun3_after  --overrepresentation_analysis --thread 16 > ${scripts}KorefMicrobiomeShotgun3_fastp_after.sh
qsub -pe smp 16 -cwd ${scripts}KorefMicrobiomeShotgun3_fastp_after.sh


END

#3. bowtie2 host removal

samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2 KorefMicrobiomeShotgun3"
output_1=${output}/1_fastp_afterTrimming/
output_2=${output}/2_bowtie_rerun/
mkdir -p $output_2

echo "start bowtie"
for sample in $samples
do
  mkdir -p ${output_2}${sample}
  touch  ${scripts}${sample}_bowtie.sh
  echo $bowtie2 -p 64 -1 ${output_1}${sample}/${sample}_1.fq.gz -2 ${output_1}${sample}/${sample}_2.fq.gz -x $hg38 -S ${output_2}${sample}/${sample}_mapped_and_unmapped.sam > ${scripts}${sample}_bowtie.sh
  echo $samtools view -bS -@ 64 ${output_2}${sample}/${sample}_mapped_and_unmapped.sam \> ${output_2}${sample}/${sample}_mapped_and_unmapped.bam >> ${scripts}${sample}_bowtie.sh
  echo $samtools view -b -f 12 -F 256 -@ 64 ${output_2}${sample}/${sample}_mapped_and_unmapped.bam \> ${output_2}${sample}/${sample}_bothEndsUnmapped.bam >> ${scripts}${sample}_bowtie.sh
  echo $samtools sort -n -@ 64 ${output_2}${sample}/${sample}_bothEndsUnmapped.bam -o ${output_2}${sample}/${sample}_bothEndsUnmapped_sorted.bam >> ${scripts}${sample}_bowtie.sh
  echo $bedtools bamtofastq -i ${output_2}${sample}/${sample}_bothEndsUnmapped_sorted.bam -fq ${output_2}${sample}/${sample}_host_removed_r1.fq -fq2 ${output_2}${sample}/${sample}_host_removed_r2.fq >> ${scripts}${sample}_bowtie.sh
  echo gzip ${output_2}${sample}/${sample}_host_removed_r1.fq >> ${scripts}${sample}_bowtie.sh
  echo gzip ${output_2}${sample}/${sample}_host_removed_r2.fq >> ${scripts}${sample}_bowtie.sh
  qsub -pe smp 64 -cwd ${scripts}${sample}_bowtie.sh

done
echo "done"

#4. quast