#!/usr/bin/env bash

# for RNAseq analysis:
conda create --name RNAseq python=3.9
# OR in a specific location
#conda create --prefix /Volumes/Seagate/opt/miniconda3/env/salmon python=3.8
# to remove
#conda remove --name RNAseq --all


CONDA_PACKAGE=`which conda`
if [ ! ${CONDA_PACKAGE} ]; then
	echo "Check if anaconda or miniconda is installed in the home directory";
    exit
fi

source ${CONDA_PACKAGE%conda}activate RNAseq
# the equaivalent will be:
# conda activate RNAseq

if [ -z ${CONDA_ACTIVATE} ]; then
  echo "export CONDA_ACTIVATE=${CONDA_PACKAGE%conda}activate" >> ~/.bashrc
fi

# OR in a specific location
conda install --update-deps cutadapt fastqc wiggletools fastp -c bioconda -c conda-forge
conda install --update-deps -c bioconda htslib samtools bcftools ucsc-wigtobigwig gffread 
#conda install --update-deps -c conda-forge ncurses

salmon_ver=1.9.0
conda install --update-deps -c bioconda salmon=$salmon_ver

STAR_ver=2.7.10b
conda install --update-deps -c bioconda star=$STAR_ver

# additional software
SOFTWARE_DIR=${HOME}/custom_software
mkdir -p $SOFTWARE_DIR
cd ${SOFTWARE_DIR}

# If conda installation doesn't work:
if ! [ -x "$(command -v salmon)" ]; then
  wget  https://github.com/COMBINE-lab/salmon/releases/download/v${salmon_ver}/salmon-${salmon_ver}_linux_x86_64.tar.gz
  tar -xzvf salmon-${salmon_ver}_linux_x86_64.tar.gz
  if [ `grep -c "salmon" ${HOME}/.bashrc` -lt 1 ]; then
    echo 'export PATH=${PATH}:'${SOFTWARE_DIR}/salmon-${salmon_ver}_linux_x86_64/bin >> ${HOME}/.bashrc
  fi
fi

# https://github.com/alexdobin/STAR
# Get latest STAR source from releases
if ! [ -x "$(command -v STAR)" ]; then
   wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_ver}.tar.gz
   tar -xzf ${STAR_ver}.tar.gz
   #cd STAR-${STAR_ver}/source
   if [ `grep -c "STAR-${STAR_ver}" ${HOME}/.bashrc` -lt 1 ];  then
     echo 'export PATH=${PATH}:'${HOME}/custom_software/STAR-${STAR_ver}/bin/Linux_x86_64_static/ >> ~/.bashrc
   fi
fi



if [ -x "$(command -v module)"  ]; then
  R_ver=4.1.0
  module load R/${R-ver};
  Rscript ./install_dependencies.R $R_ver
else
  echo "no R module"
  R_ver=4.2.1
  runrscript${R_ver}.sh ${PWD}/install_dependencies.R $R_ver
fi


