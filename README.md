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

    #fastqFile,sampleID,strain,AIDgene,TIR1,Auxin,Drug
 
If you have PE sequencing, use the following headers:
    
    #fastqFile1,fastqFile2,sampleID,strain,AIDgene,TIR1,Auxin,Drug
    
Names of the fastqFiles (with full path) and the sampleID (sampleID must be unique to each row) are essential for mapping. The other fields are used by DESeq2 to create appropriate comparison groupings and can be changed according to your data.

### 1.1 Map RNA-seq reads

Open the _01_mapRNA.sh file and change the #SBATCH --array=1-24%5 line to reflect the number of samples in your <fastqList>.csv file (here 24 files which will be processed in batches of 5 so as not to overload the server).

Make sure fastqFileList variable is set to your <fastqList>.csv file name
Make sure the genomeVer and GENOME_DIR variables are correctly set (same as in the indexing script)

In SLURM environment use the following command:

```
sbatch _01_mapRNAreads.sh 
```    

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

If you have errors with the R library installation you might have a clash of versions - you might need to remove the the libraries in the .libPaths() location and reinstall.
