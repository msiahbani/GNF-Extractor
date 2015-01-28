#!/bin/bash
#PBS -m ea
#PBS -M maryam84siabani@gmail.com
#PBS -l walltime=1:00:00
#PBS -l mem=1gb

cd /home/msiahban/git-folder/GNF-Extractor/src/training-scripts

SELL=bash
echo "source /cs/natlang-sw/`uname -s`-`uname -i`/Modules/3.2.6/init/`echo ${SELL} | awk -F "/" '{print $NF}'`"
#echo "source /cs/natlang-sw/`uname -s`-`uname -i`/Modules/3.2.6/init/`echo ${SELL} | awk -F "/" '{print $NF}'`"
source /home/msiahban/`uname -s`-`uname -i`/Modules/3.2.6/init/`echo ${SELL} | awk -F "/" '{print $NF}'`
module load use.hpc

#module load NL/LM/NGRAM-SWIG-SRILM/NGRAM
/usr/bin/env perl gnf_xtract.pl training_cn-en.config

