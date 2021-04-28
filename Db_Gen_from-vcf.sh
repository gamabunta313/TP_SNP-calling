
# $DB_V_OLD is the custom database generated from 362 TP53-WT patients. It comprises VAF observed for each TP53 positions for all the WT samples
# $DB_V_OLD is used for error modelling procedures 
# $DB_C_OLD is the custom database generated from 362 TP53-WT patients. It comprises the counts observed for each TP53 positions for all the WT samples
# $DB_C_OLD is used for the Fisher's exact test

export DIR=$PWD

scriptdir=$"/path/to/scripts"

source $scriptdir/config.txt

mkdir $DIR/DB_Generation ; DB=$DIR/DB_Generation
mkdir $DIR/RESULTS ; R=$DIR/RESULTS

for names in $DIR/V4/*_GATK-annotated.vcf
do
echo $names >> $DB/lista_vcf.list
done

cd $DB

java -Xmx8g -jar $GATK -T CombineVariants -R $REF -V lista_vcf.list -genotypeMergeOptions UNIQUIFY -o Combined_GATK.vcf

awk '$0!~/^#/' Combined_GATK.vcf | awk 'length($3)==1 && length($4)==1' | awk '{print $1,$2}' > ALL_Mutations.txt

for n in A C T G
do
awk -v x="$n" '{print $0,x}' ALL_Mutations.txt > ${n}_ALL_Mutations.txt
done

cat A_ALL_Mutations.txt C_ALL_Mutations.txt G_ALL_Mutations.txt T_ALL_Mutations.txt > ACGT_ALL_Mutations.txt

awk -F " " 'NR==FNR{a[$1,$2]=1;next}($1,$2) in a{print $0}' ACGT_ALL_Mutations.txt ${DB_C_OLD} > $R/DB_COUNT_Custom.tab
awk -F " " 'NR==FNR{a[$1,$2]=1;next}($1,$2) in a{print $0}' ACGT_ALL_Mutations.txt ${DB_V_OLD} > $R/DB_VAF_Custom.tab
cd ..
cd $DIR/SB_SNP_G

Rscript $Rscripts/Merge-table.R 
cp Good_Mutations.tab $R

cd $R

Rscript $Rscripts/Error_Modelling_Fisher-Test_TP53.R Good_Mutations.tab DB_COUNT_Custom.tab

Rscript $Rscripts/Error_Modelling_TP53_Distribution.R Good_Mutations.tab DB_VAF_Custom.tab

