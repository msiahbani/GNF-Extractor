#! /bin/tcsh 

/bin/rm -rf /global/scratch/msiahban/LR-Hiero/cn-en/Model/training/scfg-rules/devset-temp/
/bin/rm -rf /global/scratch/msiahban/LR-Hiero/cn-en/Model/training/scfg-rules/testset-temp/
/bin/rm /global/scratch/msiahban/LR-Hiero/cn-en/Model/scripts/*
unlink /global/scratch/msiahban/LR-Hiero/cn-en/Model/delete_temp_files.sh
rm -rf /local-scratch/${PBS_JOBID}
