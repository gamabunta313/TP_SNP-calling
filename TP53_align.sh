export DIR=$PWD

mkdir $DIR/BAM

B=$DIR/BAM

scriptdir=$"/dir/to/scripts"

# config file in which absolute paths are reported

source $scriptdir/config.txt

for file in $DIR/*_L001_R1_001.fastq.gz

do

R1=$(readlink -f "$file")
R2=${file%%_L001_R1_001.fastq.gz}"_L001_R2_001.fastq.gz"

filename=$(basename $file _L001_R1_001.fastq.gz)

O=$DIR/${filename}


mkdir $O ; cd $O

gunzip -c $R1 > ${O}/${filename}_R1.fastq
gunzip -c $R2 > ${O}/${filename}_R2.fastq

bwa mem -t 8 $REF ${filename}_R1.fastq ${filename}_R2.fastq | samtools sort | samtools view -Sb > ${filename}.bam

samtools index ${filename}.bam

mv ${filename}.bam* $B

cd ..
rm -rf $O

done
