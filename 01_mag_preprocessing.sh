#!/bin/bash


### paths ###

rawdata="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/raw_data/"
output="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data"
hg38="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/GRCh38_noalt_as"

samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2"
sample1="KorefMicrobiomeShotgun1"
sample2="KorefMicrobiomeShotgun2"

### programs ###

fastp="/BiO/Access/ehojune/anaconda3/bin/fastp"
megahit="/BiO/Access/ehojune/anaconda3/bin/megahit"
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py"
quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"
bowtie2="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bowtie2-2.5.1-mingw-x86_64/bowtie2"
samtools="/BiO/Access/ehojune/anaconda3/bin/samtools"


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

END
#2. fastp after setting trimming criteria

output_1=${output}/1_fastp_afterTrimming/
mkdir -p $output_1

for sample in $samples
do
  mkdir -p ${output_1}${sample}
done

echo "start fastp after"

echo $fastp -i ${rawdata}${sample1}/${sample1}_1.fq.gz -I ${rawdata}${sample1}/${sample1}_2.fq.gz -o ${output_1}${sample1}/${sample1}_1.fq.gz -O ${output_1}${sample1}/${sample1}_2.fq.gz --html ${output_1}${sample1}/${sample1}.html --json ${output_1}${sample1}/${sample1}.json --report_title ${sample1}_after  --overrepresentation_analysis --trim_front1 15 --trim_front2 15 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20


echo $fastp -i ${rawdata}${sample2}/${sample2}_1.fq.gz -I ${rawdata}${sample2}/${sample2}_2.fq.gz -o ${output_1}${sample2}/${sample2}_1.fq.gz -O ${output_1}${sample2}/${sample2}_2.fq.gz --html ${output_1}${sample2}/${sample2}.html --json ${output_1}${sample2}/${sample2}.json --report_title ${sample2}_after  --overrepresentation_analysis --trim_front1 10 --trim_front2 10 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20

mkdir -p ${output_1}KorefMicrobiomeShotgun3
cat ${output_1}${sample1}/${sample1}_1.fq.gz ${output_1}${sample2}/${sample2}_1.fq.gz > ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_1.fq.gz
cat ${output_1}${sample1}/${sample1}_2.fq.gz ${output_1}${sample2}/${sample2}_2.fq.gz > ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_2.fq.gz
$fastp -i ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_1.fq.gz -I ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_catted_2.fq.gz -o ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_1.fq.gz -O ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_2.fq.gz --html ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.html --json ${output_1}KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.json --report_title KorefMicrobiomeShotgun3_after  --overrepresentation_analysis --thread 16

echo "done"


#3. bowtie2 host removal

samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2 KorefMicrobiomeShotgun3"
output_2=${output}/2_bowtie/
mkdir -p $output_2
unzip $hg38.zip

echo "start fastp before"
for sample in $samples
do
  mkdir -p ${output_2}${sample}
  $bowtie2 -p 64 -1 ${output_1}${sample}/${sample}_1.fq.gz -2 ${output_1}${sample}/${sample}_2.fq.gz -x $hg38 --un-conc-gz ${output_2}${sample}/${sample}_mapped_and_unmapped.sam
  $samtools view -bS ${output_2}${sample}/${sample}_mapped_and_unmapped.sam > ${output_2}${sample}/${sample}_mapped_and_unmapped.bam
done
echo "done"









#4. quast