#!/bin/bash
#PBS -l arch=x86_64
#PBS -l ncpus=1
#PBS -m ea
#PBS -M maryam84siabani@gmail.com
#PBS -l walltime=24:00:00
#PBS -l mem=40gb

cd /cs/natlang-projects/LR-Hiero/Hiero-PhrXtctr/training-scripts/

SELL=bash
source /cs/natlang-sw/`uname -s`-`uname -i`/Modules/3.2.6/init/`echo ${SELL} | awk -F "/" '{print $NF}'`
module load use.hpc
module load NL/LM/NGRAM-SWIG-SRILM/NGRAM
module load NL/MT/MOSES/SVN_20101028
module load NL/MT/GIZA++/1.0.3 
perl phr_xtract.pl training.config
