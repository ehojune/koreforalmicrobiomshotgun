/BiO/Access/ehojune/anaconda3/bin/kraken2 --db /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/maxikraken2_1903_140GB --gzip-compressed --threads 50 --report /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_2_species_afterRemoval/KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1.kraken2.report.txt --output /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_2_species_afterRemoval/KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1.kraken2.tsv --paired /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_rerun//KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1_host_removed_r1.fq.gz /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_rerun//KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1_host_removed_r2.fq.gz
/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/Krona/KronaTools/bin/ktImportTaxonomy -q 2 -t 3 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_2_species_afterRemoval/KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1.kraken2.tsv -o /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/4_2_species_afterRemoval/KorefMicrobiomeShotgun1/KorefMicrobiomeShotgun1.kronaReport.html
