export PATH=/BiO/Access/ehojune/anaconda3/bin/python:$PATH
export PYTHONPATH=$PYTHONPATH:/BiO/Access/ehojune/anaconda3/bin/python
/BiO/Access/ehojune/anaconda3/bin/python /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.2-Linux/bin/metaspades.py --meta -1 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_modified/KorefMicrobiomeShotgun2/KorefMicrobiomeShotgun2_host_removed_r1.fq.gz -2 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_bowtie_modified/KorefMicrobiomeShotgun2/KorefMicrobiomeShotgun2_host_removed_r2.fq.gz -o /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/3_2_metaspades/KorefMicrobiomeShotgun2/KorefMicrobiomeShotgun2.metaspades_asm -t 100 --memory 3000
