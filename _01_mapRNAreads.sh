#! /bin/bash
#SBATCH --mail-user=jennifer.semple@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-24%5   #1-20%10

#source $HOME/.bashrc
source ${CONDA_ACTIVATE} RNAseq

echo "current date"
date
echo "current git branch"
git branch
echo "current git version"
git log -1 --format="%H"


###############################
########### VARIABLES #########
###############################

# table of fastqFiles
fastqFileList=./fastqList.csv
echo "this is $fastqFileList"
# The following command can be useful increating this file. realpath simplfies relative paths and paste puts paired end reads in columns:
# realpath ../../bisiaka/RNAseq/rawdata/X204SC21091961-Z01-F002/raw_data/*/*.fq.gz | paste -d"," - - > fastqFiles.csv
# then add the sampleID from the folder name:
# cut -d"/" -f10 fastqFiles.csv | paste -d"," fastqFiles.csv - > fastqList.csv
# rm fastqFiles.csv

genomeVer=WS285
GENOME_DIR=${HOME}/genomeData/${genomeVer}
chromSizesFile=${GENOME_DIR}/${genomeVer}.chrom.sizes

echo "GENOME_DIR=$GENOME_DIR"

mRNAindex=${GENOME_DIR}/${genomeVer}_mRNA_index
geneindex=${GENOME_DIR}/${genomeVer}_gene_index
#ncRNAindex=${GENOME_DIR}/${genomeVer}_ncRNA_index
#pseudoIndex=${GENOME_DIR}/${genomeVer}_pseudogenic_index
#tnIndex=${GENOME_DIR}/${genomeVer}_transposon_index

# check if paired end or single end
grep ",fastqFile2," $fastqFileList
if [ "$?" == 0  ]
then
    isPE=true
    echo "paired end sequences"
else
    echo "running single end commands"
fi


# output folder
WORK_DIR=${PWD}

# this selects which row of the fastqList.csv will be run
i=${SLURM_ARRAY_TASK_ID}

nThreads=${SLURM_CPUS_PER_TASK}


# get column data from file (contains all batches)
fastqFiles1=(`cut -f1 -d "," $fastqFileList`)
sampleNames=(`cut -f3 -d "," $fastqFileList`)
if [ "$isPE" == "true" ]; then
    fastqFiles2=(`cut -f2 -d "," $fastqFileList`)
else
    sampleNames=(`cut -f2 -d "," $fastqFileList`)
fi

# get line data for single batch
fastqFile1=${fastqFiles1[$i]}
if [ "$isPE" == "true" ]; then
    fastqFile2=${fastqFiles2[$i]}
fi
sampleName=${sampleNames[$i]}
echo "sampleName=$sampleName"

# not sure what this is?
#if [[ $fastqFile2 == "#*" ]]; then
#  echo "skip commented line: $fastqFile"
#  exit 0
#fi

# optional different baseName for SE and PE:
if [ "$isPE" == "true" ]; then
    baseName=${sampleName}
else
    baseName=${sampleName}
fi

echo baseName is $baseName

#####################
# mapping RNAseq data
#####################


#########################################################
#### get initial read stats                            ##
#########################################################
##
###run fastqc on sequences
#mkdir -p ${WORK_DIR}/qc/rawData
#fastqc ${fastqFile1} -t $nThreads -o ${WORK_DIR}/qc/rawData
#if [ "$isPE" == "true" ]; then
#  fastqc ${fastqFile2} -t $nThreads -o ${WORK_DIR}/qc/rawData
#fi
#
########################################################
### trim adaptors with cutadapt                       ##
########################################################
#
## use cutadapt to trim # cut Illumina adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#mkdir -p ${WORK_DIR}/cutadapt
#if [ "$isPE" == "true" ]; then
#     cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -m 15 -o ${WORK_DIR}/cutadapt/${baseName}_R1.fastq.gz -p ${WORK_DIR}/cutadapt/${baseName}_R2.fastq.gz -j $nThreads ${fastqFile1} ${fastqFile2}
#else  
#     cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 15 -o ${WORK_DIR}/cutadapt/${baseName}_R1.fastq.gz -j $nThreads ${fastqFile1}
#fi
##redo fastQC on trimmed reads
#mkdir -p ${WORK_DIR}/qc/cutadapt
#fastqc ${WORK_DIR}/cutadapt/${baseName}_R1.fastq.gz -t $nThreads -o ${WORK_DIR}/qc/cutadapt
#if [ "$isPE" == "true" ]; then
#  fastqc ${WORK_DIR}/cutadapt/${baseName}_R2.fastq.gz -t $nThreads -o ${WORK_DIR}/qc/cutadapt
#fi

######################################################
## test fastp                                       ##
######################################################
mkdir -p ${WORK_DIR}/fastp
mkdir -p ${WORK_DIR}/qc/fastp
if [ "$isPE" == "true" ]; then
  fastp -i ${fastqFile1} -I ${fastqFile2} -o ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz -O ${WORK_DIR}/fastp/${baseName}_R2.fastq.gz --thread $nThreads --detect_adapter_for_pe -j ${WORK_DIR}/qc/fastp/${baseName}_fastp.json -h ${WORK_DIR}/qc/fastp/${baseName}_fastp.html --failed_out  ${WORK_DIR}/fastp/${baseName}_failedReads.txt -R "${baseName} fastp report"
else
  fastp -i ${fastqFile1} -o ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz --thread $nThreads -j ${WORK_DIR}/qc/fastp/${baseName}_fastp.json -h ${WORK_DIR}/qc/fastp/${baseName}_fastp.html --failed_out  ${WORK_DIR}/fastp/${baseName}_failedReads.txt -R "${baseName} fastp report"
fi

multiqc -n multiqc_fastp.html -o ./qc ./qc/fastp


#######################################################
## Align to genome with STAR                         ##
#######################################################

normalisation=RPM #None or RPM
strandedness=Unstranded #Stranded or Unstranded

# align to genome
echo "aligning to genome with STAR..."
mkdir -p ${WORK_DIR}/bamSTAR

if [ "$isPE" == "true" ]; then
    STAR --genomeDir ${GENOME_DIR}  --readFilesIn ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz ${WORK_DIR}/fastp/${baseName}_R2.fastq.gz --readFilesCommand gunzip -c --outFileNamePrefix ${WORK_DIR}/bamSTAR/${baseName}_ --runThreadN $nThreads --outSAMmultNmax 1 --alignIntronMax 500 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outWigType wiggle --outWigStrand $strandedness --outWigNorm $normalisation
else
    STAR --genomeDir ${GENOME_DIR}  --readFilesIn ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz --readFilesCommand gunzip -c --outFileNamePrefix ${WORK_DIR}/bamSTAR/${baseName}_ --runThreadN $nThreads --outSAMmultNmax 1 --alignIntronMax 500 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outWigType wiggle --outWigStrand $strandedness --outWigNorm $normalisation
fi

samtools index ${WORK_DIR}/bamSTAR/${baseName}_Aligned.sortedByCoord.out.bam

wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_F_UniqueMultiple_${normalisation}.bw
wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_F_Unique_${normalisation}.bw

if [ "$strandedness" == "Stranded" ]; then
    wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_R_UniqueMultiple_${normalisation}.bw
    wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_R_Unique_${normalisation}.bw
fi

if [ -e "${WORK_DIR}/bamSTAR/${baseName}_F_UniqueMultiple_${normalisation}.bw" ]; then
    rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str1.out.wig
    rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str1.out.wig
    if [ "$strandedness" == "Stranded" ]; then
        rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str2.out.wig
        rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str2.out.wig
    fi
fi

multiqc -n multiqc_STAR.html -o ./qc ./bamSTAR


######################################################
# Count reads with Salmon                           ##
######################################################
echo "Count reads with Salmon..."


######## NOTE: I do not have estimates for --fldMean and --fldSD as i have no access to the bioanalyser files #########
#######  therefore the quantification based on the effective transcript length will be wrong!!!! Not sure how important this is?  ######

#${SALMON_SING}
if [ "$isPE" == "true" ]; then
  echo "quantify mRNA transcripts: ${WORK_DIR}/salmon/mRNA/${baseName}"
  salmon quant -i ${mRNAindex} -l A -1 ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz -2 ${WORK_DIR}/fastp/${baseName}_R2.fastq.gz  --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/mRNA/${baseName} --seqBias --gcBias --numBootstraps 100 --writeUnmappedNames
  echo "quantify genes: ${WORK_DIR}/salmon/gene/${baseName}"
  salmon quant -i ${geneindex} -l A -1 ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz -2 ${WORK_DIR}/fastp/${baseName}_R2.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/gene/${baseName} --seqBias --gcBias --numBootstraps 100 --writeUnmappedNames
else
  echo "quantify mRNA transcripts: ${WORK_DIR}/salmon/mRNA/${baseName}"
  salmon quant -i ${mRNAindex} -l A -r ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/mRNA/${baseName} --seqBias --gcBias --numBootstraps 100 --writeUnmappedNames
  echo "quantify genes: ${WORK_DIR}/salmon/gene/${baseName}"
  salmon quant -i ${geneindex} -l A -r ${WORK_DIR}/fastp/${baseName}_R1.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/gene/${baseName} --seqBias --gcBias --numBootstraps 100 --writeUnmappedNames
fi

# so far not needed
#echo "quantify ncRNA transcripts: ${WORK_DIR}/salmon/ncRNA/${baseName}"
##${SALMON_SING}
#salmon quant -i ${ncRNAindex} -l A -r ${WORK_DIR}/fastp/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/ncRNA/${baseName} --seqBias --gcBias --numBootstraps 100
#
#echo "quantify pseudoRNA transcripts: ${WORK_DIR}/salmon/pseudoRNA/${baseName}"
##${SALMON_SING}
#salmon quant -i ${pseudoIndex} -l A -r ${WORK_DIR}/fastp/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/pseudoRNA/${baseName} --seqBias --gcBias --numBootstraps 100
#
#echo "quantify TnRNA transcripts: ${WORK_DIR}/salmon/tnRNA/${baseName}"
##${SALMON_SING}
#salmon quant -i ${tnIndex} -l A -r ${WORK_DIR}/fastp/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/tnRNA/${baseName} --seqBias --gcBias --numBootstraps 100


multiqc -n multiqc_salmon-mRNA.html -o ./qc ./salmon/mRNA
multiqc -n multiqc_salmon-gene.html -o ./qc ./salmon/gene
