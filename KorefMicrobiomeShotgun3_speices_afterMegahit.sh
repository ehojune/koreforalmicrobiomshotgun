/BiO/Access/ehojune/anaconda3/bin/kraken2 --db /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/maxikraken2_1903_140GB --threads 50 --report /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_3_species_afterMegahit/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.kraken2.report.txt --output /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_3_species_afterMegahit/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.kraken2.tsv /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/3_1_megahit//KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.megahit_asm/final.contigs.fa
/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_3_species_afterMegahit/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.kraken2.tsv -o /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_3_species_afterMegahit/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.kronaReport.html