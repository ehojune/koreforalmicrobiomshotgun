export PATH=/BiO/Access/ehojune/anaconda3/bin/python:$PATH
export PYTHONPATH=$PYTHONPATH:/BiO/Access/ehojune/anaconda3/bin/python
/BiO/Access/ehojune/anaconda3/bin/python /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_rerun/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_host_removed_r1.fq.gz -2 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_rerun/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_host_removed_r2.fq.gz -o /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/3_1_megahit/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.megahit_asm --memory 0.1 -t 60
