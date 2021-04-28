
export DIR=$PWD

mkdir vcf; mkdir $DIR/vcf/V1 $DIR/vcf/V2 $DIR/vcf/V3 $DIR/vcf/V4 $DIR/vcf/TABLE

scriptdir=$"/path/to/script"

# config file in which absolute paths are reported

source $scriptdir/config.txt

for file in ${DIR}/*.bam
do
input=$(readlink -f "$file")
filename=$(basename $input .bam)
OUTDIR=${file%%.bam}

mkdir $OUTDIR ; cd $OUTDIR
cp $input $OUTDIR
samtools index $OUTDIR/${filename}.bam

$FREEBAYES -f $REF -b ${filename}.bam \
--theta 0.001 -K -e 1000 -F 0.003 -C 8 \
--min-mapping-quality 1 --min-base-quality 30 \
--use-duplicate-reads -t $BED> ${filename}.vcf

cp ${filename}.vcf $DIR/vcf/V1/

java -Xmx8g -jar $GATK \
-T SelectVariants \
-R $REF \
-L $BED \
--variant ${filename}.vcf \
-o ${filename}_GATK-filtered.vcf

sed -i 's/0\/0/0\/1/g' ${filename}_GATK-filtered.vcf

cp ${filename}_GATK-filtered.vcf $DIR/vcf/V2/

java -jar $snpEFFdir/snpEff.jar \
-c $snpEFFdir/snpEff.config \
-v -o gatk hg19 -canon ${filename}_GATK-filtered.vcf > ${filename}_snpEFF.vcf

cp ${filename}_snpEFF.vcf $DIR/vcf/V3/

java -Xmx8g -jar $GATK \
-T VariantAnnotator \
-R $REF \
-A SnpEff --variant ${filename}_GATK-filtered.vcf \
--snpEffFile ${filename}_snpEFF.vcf \
-L $BED \
-o ${filename}_GATK-annotated.vcf

$bcftools norm -m -both ${filename}_GATK-annotated.vcf | \
$vt  decompose_blocksub -  > $DIR/vcf/V4/${filename}_GATK-annotated.vcf

$bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\t%TYPE\n]' $DIR/vcf/V4/${filename}_GATK-annotated.vcf >$DIR/vcf/TABLE/${filename}_GATK-annotated.txt


cd ..
rm -rf $OUTDIR

done
