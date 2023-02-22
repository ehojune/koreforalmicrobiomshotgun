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
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.2-Linux/bin/metaspades.py"

quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"
bowtie2="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bowtie2-2.5.1/bowtie2"
samtools="/BiO/Access/ehojune/anaconda3/bin/samtools"
bedtools="/BiO/Access/ehojune/anaconda3/bin/bedtools"
python="/BiO/Access/ehojune/anaconda3/bin/python"

### commands ###








#3_1. megahit


output_2=${output}2_bowtie_modified/
output_3_1=${output}3_1_megahit/

echo "start megahit"
for sample in $samples
do
  mkdir -p ${output_3_1}${sample}
  touch  ${scripts}${sample}_megahit.sh
  echo export PATH=${python}:\$PATH > ${scripts}${sample}_megahit.sh
  echo export PYTHONPATH=\$PYTHONPATH:${python} >> ${scripts}${sample}_megahit.sh
  echo $megahit -1 ${output_2}${sample}/${sample}_host_removed_r1.fq.gz -2 ${output_2}${sample}/${sample}_host_removed_r2.fq.gz -o ${output_3_1}${sample}/${sample}.megahit_asm --memory 0.1 -t 60 >> ${scripts}${sample}_megahit.sh
  qsub -cwd -pe smp 60 ${scripts}${sample}_megahit.sh
done
echo "done"





:<<'END'

#3_2. metaspades


output_2=${output}2_bowtie_modified/
output_3_2=${output}3_2_metaspades/
mkdir ${output_3_2}

echo "start metaspades"
for sample in $samples
do
  mkdir -p ${output_3_2}${sample}
  touch  ${scripts}${sample}_metaspades.sh
  echo export PATH=${python}:\$PATH > ${scripts}${sample}_metaspades.sh
  echo export PYTHONPATH=\$PYTHONPATH:${python} >> ${scripts}${sample}_metaspades.sh
  echo $python $spades --meta -1 ${output_2}${sample}/${sample}_host_removed_r1.fq.gz -2 ${output_2}${sample}/${sample}_host_removed_r2.fq.gz -o ${output_3_2}${sample}/${sample}.metaspades_asm -t 100 --memory 3000 >> ${scripts}${sample}_metaspades.sh
  #qsub -cwd -pe smp 100 ${scripts}${sample}_metaspades.sh
done
echo "done"

END
#4. quast