#!/bin/bash


### paths ###

rawdata="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/raw_data/"
scripts="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/scripts/"
output="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data"
hg38="GRCh38_noalt_as"
#samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2"
samples="KorefMicrobiomeShotgun1 KorefMicrobiomeShotgun2 KorefMicrobiomeShotgun3"
kraken2_db="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/maxikraken2_1903_140GB"

### programs ###

fastp="/BiO/Access/ehojune/anaconda3/bin/fastp"
megahit="/BiO/Access/ehojune/anaconda3/bin/megahit"
spades="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py"
quast="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/quast/quast.py"
bowtie2="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bowtie2-2.5.1/bowtie2"
samtools="/BiO/Access/ehojune/anaconda3/bin/samtools"
bedtools="/BiO/Access/ehojune/anaconda3/bin/bedtools"
kraken2="/BiO/Access/ehojune/anaconda3/bin/kraken2"
krona="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/Krona/KronaTools/bin/ktImportTaxonomy"


### commands ###


:<<'END'
END

#4. kraken2 + krona

output_1=${output}/1_fastp_afterTrimming/
output_2=${output}/2_bowtie_rerun/
output_3_1=${output}/3_1_megahit/
output_3_2=${output}/3_2_metaspades/


output_4_1=${output}/4_1_species_beforeRemoval/
output_4_2=${output}/4_2_species_afterRemoval/
output_4_3=${output}/4_3_species_afterMegahit/
output_4_4=${output}/4_4_species_afterMetaspades/


##4-1.species before removal

mkdir -p $output_4_1

echo "start species beforeRemoval"
for sample in $samples
do
  mkdir -p ${output_4_1}${sample}
  touch ${scripts}${sample}_speices_beforeRemoval.sh
  echo $kraken2 --db $kraken2_db --gzip-compressed --threads 50 --report ${output_4_1}${sample}/${sample}.kraken2.report.txt --output ${output_4_1}${sample}/${sample}.kraken2.tsv --paired ${output_1}/${sample}/${sample}_1.fq.gz ${output_1}/${sample}/${sample}_2.fq.gz > ${scripts}${sample}_speices_beforeRemoval.sh
  echo $krona -q 2 -t 3 ${output_4_1}${sample}/${sample}.kraken2.tsv -o ${output_4_1}${sample}/${sample}.kronaReport.html >> ${scripts}${sample}_speices_beforeRemoval.sh
  qsub -cwd -pe smp 100 ${scripts}${sample}_speices_beforeRemoval.sh  
done
echo "done"



mkdir -p $output_4_2

echo "start species afterRemoval"
for sample in $samples
do
  mkdir -p ${output_4_2}${sample}
  touch ${scripts}${sample}_speices_afterRemoval.sh
  echo $kraken2 --db $kraken2_db --gzip-compressed --threads 50 --report ${output_4_2}${sample}/${sample}.kraken2.report.txt --output ${output_4_2}${sample}/${sample}.kraken2.tsv --paired ${output_2}/${sample}/${sample}_host_removed_r1.fq.gz  ${output_2}/${sample}/${sample}_host_removed_r2.fq.gz > ${scripts}${sample}_speices_afterRemoval.sh
  echo $krona -q 2 -t 3 ${output_4_2}${sample}/${sample}.kraken2.tsv -o ${output_4_2}${sample}/${sample}.kronaReport.html >> ${scripts}${sample}_speices_afterRemoval.sh
  qsub -cwd -pe smp 100 ${scripts}${sample}_speices_afterRemoval.sh  
done
echo "done"




mkdir -p $output_4_3

echo "start species afterMegahit"
for sample in $samples
do
  mkdir -p ${output_4_3}${sample}
  touch ${scripts}${sample}_speices_afterMegahit.sh
  echo $kraken2 --db $kraken2_db --threads 50 --report ${output_4_3}${sample}/${sample}.kraken2.report.txt --output ${output_4_3}${sample}/${sample}.kraken2.tsv ${output_3_1}/${sample}/${sample}.megahit_asm/final.contigs.fa > ${scripts}${sample}_speices_afterMegahit.sh
  echo $krona -q 2 -t 3 ${output_4_3}${sample}/${sample}.kraken2.tsv -o ${output_4_3}${sample}/${sample}.kronaReport.html >> ${scripts}${sample}_speices_afterMegahit.sh
  qsub -cwd -pe smp 100 ${scripts}${sample}_speices_afterMegahit.sh
done
echo "done"



mkdir -p $output_4_4

echo "start species afterMetaspades"
for sample in $samples
do
  mkdir -p ${output_4_4}${sample}
  touch ${scripts}${sample}_speices_afterMetaspades.sh
  echo $kraken2 --db $kraken2_db --threads 50 --report ${output_4_4}${sample}/${sample}.kraken2.report.txt --output ${output_4_4}${sample}/${sample}.kraken2.tsv ${output_3_2}/${sample}/${sample}.metaspades_asm/scaffolds.fasta > ${scripts}${sample}_speices_afterMetaspades.sh
  echo $krona -q 2 -t 3 ${output_4_4}${sample}/${sample}.kraken2.tsv -o ${output_4_4}${sample}/${sample}.kronaReport.html >> ${scripts}${sample}_speices_afterMetaspades.sh
  qsub -cwd -pe smp 100 ${scripts}${sample}_speices_afterMetaspades.sh
done
echo "done"
