#!/bin/bash


### paths ###

rawdata="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/raw_data/"
output="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/"
scripts="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/scripts/"


samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2 KorefMicrobiomeShotgun3"
sample1="KorefMicrobiomeShotgun1"
sample2="KorefMicrobiomeShotgun2"
sample3="KorefMicrobiomeShotgun3"

### programs ###

fastp="/BiO/Access/ehojune/anaconda3/bin/fastp"
megahit="/BiO/Access/ehojune/anaconda3/bin/megahit"
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py"
quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"


### commands ###


:<<'END'

# 1. fastp before setting trimming criteria

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

echo $fastp -i ${rawdata}${sample1}/${sample1}_1.fq.gz -I ${rawdata}${sample1}/${sample1}_2.fq.gz -o ${output_1}${sample1}/${sample1}_1.fq.gz -O ${output_1}${sample1}/${sample1}_2.fq.gz --html ${output_1}${sample1}/${sample1}.html --json ${output_1}${sample1}/${sample1}.json --report_title ${sample1}_after  --overrepresentation_analysis --trim_front1 15 --trim_front2 15 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20


echo $fastp -i ${rawdata}${sample2}/${sample2}_1.fq.gz -I ${rawdata}${sample2}/${sample2}_2.fq.gz -o ${output_1}${sample2}/${sample2}_1.fq.gz -O ${output_1}${sample2}/${sample2}_2.fq.gz --html ${output_1}${sample2}/${sample2}.html --json ${output_1}${sample2}/${sample2}.json --report_title ${sample2}_after  --overrepresentation_analysis --trim_front1 8 --trim_front2 8 --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 0 --thread 16 --cut_front --average_qual 20

echo "done"







#3. megahit

output_1=${output}/1_fastp_afterTrimming/
output_2=${output}/2_1_megahit/
mkdir -p $output_2

echo "start megahit"
for sample in $samples
do
  mkdir -p ${output_2}${sample}
  touch /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/${sample}.megahit.sh
  echo $megahit -1 ${output_1}${sample}/${sample}_1.fq.gz -2 ${output_1}${sample}/${sample}_2.fq.gz -o ${output_2}${sample}/${sample}.megahit_asm --memory 0.1 -t 20> /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/${sample}.megahit.sh
done
echo "done"
END


#4. metaspades
output_1=${output}/1_fastp_afterTrimming/
output_3=${output}/2_2_metaspades/
mkdir -p $output_3

echo "start megaspades"
for sample in $samples
do
  mkdir -p ${output_3}${sample}
  touch ${scripts}/${sample}.metaspades.sh
  echo $spades --meta -1 ${output_1}${sample}/${sample}_1.fq.gz -2 ${output_1}${sample}/${sample}_2.fq.gz -o ${output_3}${sample}/${sample}.metaspades_asm -t 20 > /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/${sample}.metaspades.sh
done
echo "done"




END







#5. quast
output_2=${output}/2_1_megahit/
output_3=${output}/2_2_metaspades/
output_5=${output}/3_1_quast/

mkdir -p $output_5

echo "start quast"
for sample in $samples
do
  mkdir -p ${output_5}/megahit/${sample}
  touch ${scripts}/${sample}.megahit_quast.sh
  echo $quast ${output_2}/${sample}/${sample}.megahit_asm/final.contigs.fa -o ${output_5}/megahit/${sample} > ${scripts}/${sample}.megahit_quast.sh
  qsub -cwd -pe smp 16 ${scripts}/${sample}.megahit_quast.sh

  mkdir -p ${output_5}/metaspades/${sample}
  touch ${scripts}/${sample}.metaspades_quast.sh
  echo $quast ${output_3}/${sample}/${sample}.metaspades_asm/final.contigs.fa -o ${output_5}/megahit/${sample} > ${scripts}/${sample}.megahit_quast.sh
  qsub -cwd -pe smp 16 ${scripts}/${sample}.megahit_quast.sh
  
  
done











echo "done"

#5. quast
output_3=${output}/2_2_metaspades/
output_2=${output}/2_1_megahit/

mkdir -p $output_3

echo "start megaspades"
for sample in $samples
do
  mkdir -p ${output_3}/${sample}
  touch /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/${sample}.metaspades.sh
  echo $spades --meta -1 ${output_1}/${sample}/${sample}_1.fq.gz -2 ${output_1}/${sample}/${sample}_2.fq.gz -o ${output_3}/${sample}/${sample}.metaspades_asm -t 20 > /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/${sample}.metaspades.sh
done
echo "done"



























