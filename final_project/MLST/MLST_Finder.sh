#!/bin/sh

#****************************************************************
#  MLST_finder.2.sh
#  
#  Created by Meng Ho on 4/6/19.
#  This script is intended for the MCB432 Course as an
#  In-class shell-scripting exercise
#  Contact mho1@illinois.edu for comments.
#
# Input file is a fasta file containing multiple sequences
# Usage: MLST_finder.2.sh file.fasta
# Download database from https://pubmlst.org/data/
# For E coli 7 gene MLST profiles, and for all seven gene alleles
#****************************************************************

# Rename and makeblastdb
# Skip this line if the data blastdb has already been formatted
#for p in adk fumC gyrB icd mdh purA recA; do mv $p.tfa.txt $p.tfa; makeblastdb -in $p.tfa -dbtype nucl;done
for p in adk fumC gyrB icd mdh purA recA; do makeblastdb -in $p.tfa -dbtype nucl;done

# Run for loop from the list "adk fumC gyrB icd mdh purA recA"
# pick the best hit for assignment of the 7-gene allele
## for p in adk fumC gyrB icd mdh purA recA; do blastn -db $p.tfa -query $1 -out $p.blastout.txt -num_descriptions 5 -num_alignments 1; awk -v nam=$p '/Query=/{Qu=$2}/>/{print Qu,substr($1,2)}/No hits found/{print Qu,nam"_0"}' <$p.blastout.txt>$p.hit.txt;done
## for p in adk fumC gyrB icd mdh purA recA; do blastn -db $p.tfa -query $1 -out $p.blastout.txt -num_descriptions 5 -num_alignments 1; awk -v nam=$p '/Query=/{Qu=$2}/>/{print Qu=$2}/No hits found/{print Qu,nam"_0"}' <$p.blastout.txt>$p.hit.txt;done

for p in adk fumC gyrB icd mdh purA recA; do blastn -db $p.tfa -query $1 -out $p.blastout.txt -num_descriptions 5 -num_alignments 1; awk -v nam=$p '/Query=/{Qu=$2}/>/{print Qu,$2}/No hits found/{print Qu,nam"_0"}' <$p.blastout.txt>$p.hit.txt;done

# Step#1: Combine all *.hit.txt lists into one single list
for p in adk fumC gyrB icd mdh purA recA; do cat $p.hit.txt; done >tmp.1.txt
# Step#2: Sort the combined list
sort tmp.1.txt >tmp.2.txt

# Combine Steps #1, #2
# for p in adk fumC gyrB icd mdh purA recA; do cat $p.hit.txt; done |sort >tmp.2.txt

# Step #3: Consolidate 7-gene MLST allele asignment into single lines
awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.2.txt >tmp.3.txt

# Step #4: tabulate the digitized allelic contents
awk 'BEGIN{FS="[ ]*[_]*"}{print $1,$3,$5,$7,$9,$11,$13,$15}' <tmp.3.txt >tmp.4.txt

# Combine Steps #3, #4
# awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.2.txt| awk 'BEGIN{FS="[ ]*[_]*"}{print $1,$3,$5,$7,$9,$11,$13,$15}' >tmp.4.txt

# Step#5: Combine Steps #3, #4 oupput as tab-delimited text
awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.2.txt| awk 'BEGIN{FS="[ ]*[_]*"}{print $1"\t"$3"\t"$5"\t"$7"\t"$9"\t"$11"\t"$13"\t"$15}'>tmp.5.txt

# Step#6: Check for ST type according to allelic asignment
awk 'BEGIN{while((getline<"ecoli.txt")>0){n++;st[n]=$1; adk[n]=$2; fumC[n]=$3; gyrB[n]=$4; icd[n]=$5; mdh[n]=$6; purA[n]=$7; recA[n]=$8}; print "Accession\tst\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA"}{for(i=2;i<=n;i++){if($2==adk[i] && $3==fumC[i] && $4==gyrB[i] && $5==icd[i] && $6==mdh[i] && $7==purA[i] && $8==recA[i]){print $1"\t"st[i]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;break}}}' <tmp.5.txt >tmp.6.txt

# Combine Steps #5, #6
#awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.2.txt| awk 'BEGIN{FS="[ ]*[_]*"}{print $1"\t"$3"\t"$5"\t"$7"\t"$9"\t"$11"\t"$13"\t"$15}'|awk 'BEGIN{while((getline<"ecoli.txt")>0){n++;st[n]=$1;adk[n]=$2;fumC[n]=$3;gyrB[n]=$4;icd[n]=$5;mdh[n]=$6;purA[n]=$7;recA[n]=$8};print "Accession\tst\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA"}{for(i=2;i<=n;i++){if($2==adk[i] && $3==fumC[i] && $4==gyrB[i] && $5==icd[i] && $6==mdh[i] && $7==purA[i] && $8==recA[i]){print $1"\t"st[i]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;break}}}'  >tmp.6.txt

# Combine Steps #1 through #6
#for p in adk fumC gyrB icd mdh purA recA; do cat $p.hit.txt; done |sort | awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'|awk 'BEGIN{FS="[ ]*[_]*"}{print $1"\t"$3"\t"$5"\t"$7"\t"$9"\t"$11"\t"$13"\t"$15}'|awk 'BEGIN{while((getline<"ecoli.txt")>0){n++;st[n]=$1;adk[n]=$2;fumC[n]=$3;gyrB[n]=$4;icd[n]=$5;mdh[n]=$6;purA[n]=$7;recA[n]=$8};print "Accession\tst\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA"}{for(i=2;i<=n;i++){if($2==adk[i] && $3==fumC[i] && $4==gyrB[i] && $5==icd[i] && $6==mdh[i] && $7==purA[i] && $8==recA[i]){print $1"\t"st[i]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;break}}}'>tmp.6.txt

# Final MLST output
rm final.txt; awk '{n++}n==1{print >"final.txt"}n>1{print}'<tmp.6.txt |sort -k2n >>final.txt

echo
cat final.txt
mv final.txt $1.result.txt
rm tmp*

