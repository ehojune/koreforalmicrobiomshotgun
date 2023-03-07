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
python="/BiO/Access/ehojune/anaconda3/bin/python"
busco="/BiO/Access/ehojune/anaconda3/envs/buscogogo/bin/busco"
bbwrap="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bbmap/bbwrap.sh"
bbmappileup="/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/bbmap/pileup.sh"
checkm="/BiO/Access/ehojune/anaconda3/bin/checkm"

### commands ###


:<<'END'
END



output_1=${output}/1_fastp_afterTrimming/
output_2=${output}/2_bowtie_rerun/
output_3_1=${output}/3_1_megahit/
output_3_2=${output}/3_2_metaspades/
output_4_1=${output}/4_1_species_beforeRemoval/
output_4_2=${output}/4_2_species_afterRemoval/
output_4_3=${output}/4_3_species_afterMegahit/
output_4_4=${output}/4_4_species_afterMetaspades/

output_5_1=${output}/5_1_quast_megahit/
output_5_2=${output}/5_2_quast_metaspades/

output_6_1=${output}/6_1_busco_megahit/
output_6_2=${output}/6_2_busco_metaspades/

output_7_1=${output}/7_1_bbmap_megahit/
output_7_2=${output}/7_2_bbmap_metaspades/

output_8_1=${output}/8_1_checkM_megahit/
output_8_1=${output}/8_1_checkM_metaspades/

:<<'END'

#5-1. quast for megahit
mkdir -p $output_5_1

echo "start quast for megahit"
for sample in $samples
do
  mkdir -p ${output_5_1}${sample}
  touch ${scripts}${sample}_quast_megahit.sh
  echo $python $quast ${output_3_1}/${sample}/${sample}.megahit_asm/final.contigs.fa -o ${output_5_1}/${sample}/${sample} -t 16 > ${scripts}${sample}_quast_megahit.sh
  qsub -cwd -pe smp 16 ${scripts}${sample}_quast_megahit.sh
done
echo "done"



#5-2. quast for metaspades
mkdir -p $output_5_2

echo "start quast for metaspades"
for sample in $samples
do
  mkdir -p ${output_5_2}${sample}
  touch ${scripts}${sample}_quast_metaspades.sh
  echo $python $quast ${output_3_2}/${sample}/${sample}.metaspades_asm/scaffolds.fasta -o ${output_5_2}/${sample}/${sample} -t 16 > ${scripts}${sample}_quast_metaspades.sh
  qsub -cwd -pe smp 16 ${scripts}${sample}_quast_metaspades.sh
done
echo "done"









#6-1. busco for megahit
## This should be done in seperate env, manually (not using qsub but "sh [buscofile.sh]" )
## See commands under "Alternatively you can create a new environment with BUSCO installed" sentence in   https://busco.ezlab.org/busco_userguide.html 

mkdir -p $output_6_1

echo "start writing busco scripts for megahit"
for sample in $samples
do
  mkdir -p ${output_6_1}${sample}
  touch ${scripts}${sample}_busco_megahit.sh
  echo $busco -i ${output_3_1}/${sample}/${sample}.megahit_asm/final.contigs.fa -o ${output_6_1}/${sample} --auto-lineage -m genome -c 100 > ${scripts}${sample}_busco_megahit.sh

done
echo "done"


#6-2. busco for metaspades
## This should be done in seperate env, manually (not using qsub but "sh [buscoshellfile.sh]" )
## See commands under "Alternatively you can create a new environment with BUSCO installed" sentence in   https://busco.ezlab.org/busco_userguide.html 

mkdir -p $output_6_2

echo "start writing busco scripts for metaspades"
for sample in $samples
do
  mkdir -p ${output_6_2}${sample}
  touch ${scripts}${sample}_busco_metaspades.sh
  echo $busco -i ${output_3_2}/${sample}/${sample}.metaspades_asm/scaffolds.fasta -o ${output_6_2}/${sample} --auto-lineage -m genome -c 100 > ${scripts}${sample}_busco_metaspades.sh

done
echo "done"

END



#7-1. bbmap for megahit
#CheckM was installed to ehojune's conda env on shrimp, referring https://github.com/Ecogenomics/CheckM/wiki/Installation#installation-through-pip



mkdir -p $output_7_1

echo "start writing bbmap scripts for megahit"
for sample in $samples
do
  mkdir -p ${output_7_1}${sample}
  touch ${scripts}${sample}_bbmap_megahit.sh
  echo $bbwrap ref=${output_3_1}/${sample}/${sample}.megahit_asm/final.contigs.fa in=${output_2}/${sample}/${sample}_host_removed_r1.fq.gz in2=${output_2}/${sample}/${sample}_host_removed_r2.fq.gz out=${output_7_1}/${sample}/${sample}.bbmap.aln.sam.gz kfilter=22 subfilter=15 maxindel=80 >> ${scripts}${sample}_bbmap_megahit.sh
  echo $bbmappileup in=${output_7_1}/${sample}/${sample}.bbmap.aln.sam.gz out=${output_7_1}/${sample}/${sample}.cov.txt >> ${scripts}${sample}_bbmap_megahit.sh
  echo $samtools view -u -f4 ${output_7_1}/${sample}/${sample}.bbmap.aln.sam.gz | $samtools bam2fq    >> ${scripts}${sample}_bbmap_megahit.sh
  qsub -cwd -pe smp 16 ${scripts}${sample}_bbmap_megahit.sh

done
echo "done"




$bbwrap ref=SRR341725.megahit_asm/final.contigs.fa in=SRR341725_1.fastq.gz in2=SRR341725_2.fastq.gz out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80
$bbmappileup in=aln.sam.gz out=cov.txt
samtools-1.2/samtools view -u -f4 aln.sam.gz | samtools-1.2/samtools bam2fq -s unmapped.se.fq - > unmapped.pe.fq
# trying to change: use gatk instead of samtools
java -jar picard.jar qsub I=merged.dedup.realn.bam FASTQ=output.samtools.bam2fq_R1.fastq --SECOND_END_FASTQ=outout.samtools.bam2fq_R2.fastq




#7-2. bbmap for metaspades






#8-1. checkM for megahit
#CheckM was installed to ehojune's conda env on shrimp, referring https://github.com/Ecogenomics/CheckM/wiki/Installation#installation-through-pip





#8-2. checkM for metaspades



