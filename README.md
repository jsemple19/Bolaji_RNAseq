# Bolaji_RNAseq

Based on pipeline created by Jenny Semple (SMC_RNAseq) and Todor Gitchev(CeFTALL).

## Installation and preparation

This only needs to be done once

### 0.1 Clone project:

From command line:

    git clone git@github.com:CellFateNucOrg/Bolaji_RNAseq.git

    # how to update the remote url when iusing access token:
    git remote set-url origin git@github.com:CellFateNucOrg/Bolaji_RNAseq.git
    # how to see what is current set:
    git config --get remote.origin.url
    # after propperly set accessible remote url, you should be able to do: git pull 

### 0.2 Python conda:

Python Miniconda installation if needed:  
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh -b
    
    # to use conda with bash shell 
    conda init
```

### 0.3. Install environment:

Installation will be done primarily with conda you need to make sure the base environment is up to date. 

run script: 
```
    ./install_ceFTALL_env.sh
```

For R installation there is some variation between bioinformatics and ubelix clusters where you can use module load and Pertz cluster where you use a singularity image of R.
  

### 0.4. Downloads reference data and index genome/transcripts

Change location of where to put reference genome in _00_downloadAndIndexGenome.sh
run script:

    ./_00_downloadAndIndexGenome.sh

   
## Run RNA-seq analysis

### 1.1 Prepare your data

Place all desired .fastq.gz files in a folder and create a <fastqList>.csv with the column structure:

    #fastqFile,lane,sampleID,strain,AIDgene,TIR1,Auxin,Drug

and as row properly list your files or interest. Normally this is a selection from the file sampleList.csv


### 1.1 Map RNA-seq reads
Go into the src file to run the scripts.

```
cd src
```

In SLURM environment use the following command:

    sbatch --array=1-<number of fastq files>%<max number of simultaneously running jobs> sbatch_01_mapRNAreads.sh <fastqList.csv>

Examples run:
    
    # all
    sbatch  --array=1-18%9 sbatch_01_mapRNAreads.sh ./data/fastqList_t-all_ringo_kingston.csv
    # single one or selected
    sbatch  --array=7,8 sbatch_01_mapRNAreads.sh ./data/fastqList_t-all_ringo_kingston.csv

In a Unix shell console (Mac/Linux) without SLURM you can use the following command:

    <oath>/run_01_mapRNAreads.sh <fastqList.csv> <outpuFolder>

It will run for all fastq files in the list (<fastqList.csv>)

Example:

    ./src/run_01_mapRNAreads.sh ./data/ceFTALL_fastqList.csv ./out

Optionally for a single fastq file, you can use:
    
    <oath>/_01_mapRNAreads.sh <fastqList.csv> <outpuFolder> <fastq number==line number -1>
    

### 1.1 Differential Expression (DESeq) analysis



In SLURM environment, start differential expression analysis with the command: 

    sbatch <path>/sbatch_02_ceDESeqAnalysis_salmon.sh <fastqList.csv> <outpuFolder>

Example usage:

    sbatch sbatch_02_ceDESeqAnalysis_salmon.sh ./data/fastqList_t-all_ringo_kingston.csv
    sbatch sbatch_02_ceDESeqAnalysis_salmon.sh ./data/fastqList_t-all_ringo_kingston_284_285.csv

### A. Troubleshooting

Check log files and if you notice error messages with missing software reinstall the missing package. Example:
Error message:
    wiggletools: not found
    # install it in our conda environment:
    conda activate ceftall
    conda install wiggletools
    

Sometimes it's just enough to activate manually the environment before running a sbatch script:
    conda activate ceftall

