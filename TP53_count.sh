
## this script is designed to take as argument .bam files and generate a count table.
## Counting process is based on samtools 1.18 and pileup2baseindel with slight modifications
## To use the script you shoud run it in the folder in which .bam files are contained 

export DIR=$PWD

mkdir $DIR/vcf/COUNTS
mkdir $DIR/vcf/RAW_COUNTS

scriptdir=$"/path/to/scripts"

source $scriptdir/config.txt

for file in ${DIR}/*.bam
do

input=$(readlink -f "$file")
filename=$(basename $input .bam)
OUTDIR=${file%%.bam}

mkdir $OUTDIR ; cd $OUTDIR
cp $input $OUTDIR; cp ${input}.bai $OUTDIR

$samtools view -H $input | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | \
xargs -I {} -n 1 -P 8 sh -c "$samtools mpileup -ABQ0 -d 10000000000 -f $REF -r {} -l $BED $input > ${filename}_tmp.{}.mpileup"

for file in ./*${filename}*.mpileup;
do
file_name=$(basename $file .mpileup)
row=$(cat $file | wc -l)
if [ "$row" == "0" ];then
rm -rf $file
else 
awk '{gsub(/c/, "C", $3);gsub(/g/, "G", $3);gsub(/a/, "A", $3);gsub(/t/, "T", $3)}1' $file > ${file_name}.2.mpileup
fi
done
$samtools view -H $input | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | \
xargs -I {} -n 1 -P 8 sh -c "path/to/pileup2baseindel.pl -i ${filename}_tmp.{}.2.mpileup -bq 0 -prefix ${filename}_pileup2base_Q0.{}.txt"
$samtools view -H $input | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | \
xargs -I {} -n 1 -P 8 sh -c "path/to/pileup2baseindel.pl -i ${filename}_tmp.{}.2.mpileup -bq 30 -prefix ${filename}_pileup2base_Q30.{}.txt"
for fl in ./${filename}_pileup2base_Q0*
do
sed -i '1d' $fl
done
for fl in ./${filename}_pileup2base_Q30*
do
sed -i '1d' $fl
done
cat ${filename}_pileup2base_Q0* | awk -F "\t" '{OFS=FS; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' >${filename}_pileup2base_Q0.ALL.txt
cat ${filename}_pileup2base_Q30* | awk -F "\t" '{OFS=FS; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' >${filename}_pileup2base_Q30.ALL.txt

awk -F "\t" 'NR==FNR{a[$1"_"$2]=$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11;next}($1"_"$2) in a {print $0"\t"a[$1"_"$2]}' \
${filename}_pileup2base_Q30.ALL.txt ${filename}_pileup2base_Q0.ALL.txt | awk -F " " 'NR>0 ; BEGIN{OFS="\t"; print "chr","pos","REF","A0","T0","C0","G0","a0","t0","c0","g0","A30","T30","C30","G30","a30","t30","c30","g30"}' > ${filename}_pileup2base_Q0-Q30.ALL.txt

awk -F "\t" 'NR>1{
sum0=($4+$5+$6+$7+$8+$9+$10+$11)
sum30=($12+$13+$14+$15+$16+$17+$18+$19)
div1= $4 ? ($12/$4) : $12
div2= $5 ? ($13/$5) : $13
div3= $6 ? ($14/$6) : $14
div4= $7 ? ($15/$7) : $15
div5= $8 ? ($16/$8) : $16
div6= $9 ? ($17/$9) : $17
div7= $10 ? ($18/$10) : $18
div8= $11 ? ($19/$11) : $19
div9= sum0 ? (sum30/sum0) : sum30
A= ($12+$16)
T= ($13+$17)
C= ($14+$18)
G= ($15+$19)
A30= sum30 ? (A/sum30) : 0
T30= sum30 ? (T/sum30) : 0
C30= sum30 ? (C/sum30) : 0
G30= sum30 ? (G/sum30) : 0
print $1,$2,$3,div1,div2,div3,div4,div5,div6,div7,div8,div9,A,T,C,G,A30,T30,C30,G30}' ${filename}_pileup2base_Q0-Q30.ALL.txt | awk 'BEGIN{print "chr pos REF A30/A0 T30/T0 C30/C0 G30/G0 a30/a0 t30/t0 c30/c0 g30/g0 tot30/tot0 A T C G VAF-A VAF-T VAF-C VAF-G"}1' | \
awk -F " " '{gsub(" ","\t");print $0}' > ${filename}_pileup2base_Q0-Q30.ALL_Percentages.txt
cp ${filename}_pileup2base_Q0-Q30.ALL_Percentages.txt $DIR/vcf/COUNTS
cp ${filename}_pileup2base_Q30.ALL.txt $DIR/vcf/RAW_COUNTS

cd ..
rm -rf $OUTDIR
done
