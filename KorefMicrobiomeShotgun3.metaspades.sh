export PATH=/BiO/Access/ehojune/anaconda3/bin/python:$PATH
export PYTHONPATH="${PYTHONPATH}:/BiO/Access/ehojune/anaconda3/bin/python"


/BiO/Access/ehojune/anaconda3/bin/python /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/tools/SPAdes-3.15.5-Linux/bin/metaspades.py --meta -1 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/1_fastp_afterTrimming/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_1.fq.gz -2 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/1_fastp_afterTrimming/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3_2.fq.gz -o /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/ehojune/processed_data/2_2_megaspades/KorefMicrobiomeShotgun3/KorefMicrobiomeShotgun3.metaspades_asm -t 40 -m 1000
