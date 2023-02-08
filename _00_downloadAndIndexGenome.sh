#! /bin/bash
#SBATCH --mail-user=jennifer.semple@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="index_gentrome"
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G


# script to obtiain txptome fasta files from wormbase and index them for salmon
genomeVer=WS285
genomeDir=${HOME}/genomeData/${genomeVer}
genomeFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
nThreads=$SLURM_CPUS_PER_TASK

mkdir -p ${genomeDir}

source $CONDA_ACTIVATE RNAseq

#############
## index for STAR
#############
# genome sequence (only necessary for STAR)
if [ ! -f "${genomeFile}" ]; then 
  wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.${genomeVer}.genomic.fa.gz
  gunzip ${genomeFile}
fi

# create gtf annotation file using wormbase's gff file
gffFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
if [ ! -f "${gtfFile}" ]; then
  wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3.gz 
  gunzip ${gffFile}.gz
  
  # use gffread to convert gff to gtf
  b=(`basename -s .gff3 ${gffFile}`)
  gffread $gffFile -T -o ${gtfFile}
  rm $gffFile

  # need to remove wierd exons: (could maybe do this with gffread options)
  grep WormBase ${gtfFile} > ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf
  mv  ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf ${gtfFile}
fi


#index genome
#echo "indexing genome..."
STAR --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${gtfFile} --runThreadN $nThreads --genomeSAindexNbases 12 

## get chrom.sizes for wigToBigWig
if [ -x "$(command -v module)"  ]; then
  R_ver=4.1.0
  module load R/${R-ver};
  Rscript ${PWD}/getChromSizesFile.R $genomeVer $genomeDir
else
  echo "no R module"
  R_ver=4.2.1
  runrscript${R_ver}.sh ${PWD}/getChromSizesFile.R $genomeVer $genomeDir 
fi


##############
## index for Salmon
##############

# get masked genome seqence to serve as decoy
wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz
  salmon quant -i ${geneindex} -l A -r ${WORK_DIR}/cutadapt/${baseName}_R1.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/gene/${baseName} --seqBias --gcBias --numBootstraps 100

grep "^>" <(gunzip -c  ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz) | cut -d " " -f 1 > ${genomeDir}/decoys.txt
sed -i.bak -e 's/>//g' ${genomeDir}/decoys.txt

# transcript sequences (for salmon)
wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz

cat  ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${genomeDir}/mRNA_gentrome.fa.gz

salmon index -t ${genomeDir}/mRNA_gentrome.fa.gz -d ${genomeDir}/decoys.txt -p $nThreads -i ${genomeDir}/${genomeVer}_mRNA_index

rm ${genomeDir}/mRNA_gentrome.fa.gz


##wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.CDS_transcripts.fa.gz
##wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.ncRNA_transcripts.fa.gz
##wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.pseudogenic_transcripts.fa.gz
##wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.transposon_transcripts.fa.gz
#
#
##salmon index -t ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.ncRNA_transcripts.fa.gz -i ${genomeDir}/${genomeVer}_ncRNA_index
##salmon index -t ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.pseudogenic_transcripts.fa.gz -i ${genomeDir}/${genomeVer}_pseudogenic_index
##salmon index -t ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.transposon_transcripts.fa.gz -i ${genomeDir}/${genomeVer}_transposon_index


# geneIDs
wget -P ${genomeDir} ftp://ftp.wormbase.org/pub/wormbase/releases/${genomeVer}/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.${genomeVer}.geneIDs.txt.gz


# Genes for dExon/dIntron analysis:
if [ -x "$(command -v module)"  ]; then
  R_ver=4.1.0
  module load R/${R-ver};
  Rscript ${PWD}/getGeneSeqs.R $genomeVer $genomeDir
else
  echo "no R module"
  R_ver=4.2.1
  runrscript${R_ver}.sh ${PWD}/getGeneSeqs.R $genomeVer $genomeDir
fi


cat  ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genes.fa.gz ${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${genomeDir}/gene_gentrome.fa.gz

salmon index -t ${genomeDir}/gene_gentrome.fa.gz -d ${genomeDir}/decoys.txt -p $nThreads -i ${genomeDir}/${genomeVer}_gene_index

rm ${genomeDir}/gene_gentrome.fa.gz
