#!/bin/bash


### paths ###

rawdata="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/raw_data/"
scripts="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/scripts/"
output="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/"
hg38="GRCh38_noalt_as"

samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2  KorefMicrobiomeShotgun3"
sample1="KorefMicrobiomeShotgun1"
sample2="KorefMicrobiomeShotgun2"
sample3="KorefMicrobiomeShotgun3"

### programs ###

fastp="/BiO/Access/ehojune/anaconda3/bin/fastp"
megahit="/BiO/Access/ehojune/anaconda3/bin/megahit"
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py"
quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"
bowtie2="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bowtie2-2.5.1/bowtie2"
samtools="/BiO/Access/ehojune/anaconda3/bin/samtools"
bedtools="/BiO/Access/ehojune/anaconda3/bin/bedtools"
python="/BiO/Access/ehojune/anaconda3/bin/python"



output_2=${output}2_bowtie_modified/
echo "start samtools and bedtools"

for sample in $samples
do
  touch  ${scripts}${sample}_bowtie_then_modify.sh
  echo $samtools view -bS ${output_2}${sample}/${sample}_mapped_and_unmapped.sam \> ${output_2}${sample}/${sample}_mapped_and_unmapped.bam > ${scripts}${sample}_bowtie_then_modify.sh
  echo $samtools view -b -f 12 -F 256 ${output_2}${sample}/${sample}_mapped_and_unmapped.bam \> ${output_2}${sample}/${sample}_bothEndsUnmapped.bam >> ${scripts}${sample}_bowtie_then_modify.sh
  echo $samtools sort -n ${output_2}${sample}/${sample}_bothEndsUnmapped.bam ${output_2}${sample}/${sample}_bothEndsUnmapped.sorted.bam >> ${scripts}${sample}_bowtie_then_modify.sh
  echo $bedtools bamtofastq -i ${output_2}${sample}/${sample}_bothEndsUnmapped.sorted.bam -fq ${output_2}${sample}/${sample}_host_removed_r1.fq -fq2 ${output_2}${sample}/${sample}_host_removed_r2.fq >> ${scripts}${sample}_bowtie_then_modify.sh
  echo gzip ${output_2}${sample}/${sample}_host_removed_r1.fq >> ${scripts}${sample}_bowtie_then_modify.sh
  echo gzip ${output_2}${sample}/${sample}_host_removed_r2.fq >> ${scripts}${sample}_bowtie_then_modify.sh

  qsub -cwd  ${scripts}${sample}_bowtie_then_modify.sh
done
echo "done"

