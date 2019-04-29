#!/bin/sh

#****************************************************************
#  ResFinder.2.sh
#
#  Created by Meng Ho on 4/16/19.
#  This script is intended for the MCB432 Course as an
#  In-class shell-scripting exercise
#  Contact mho1@illinois.edu for comments.
#
# Input file is a fasta file containing multiple sequences
# Usage: ResFinder.2.sh file.fasta
# Download database from https://cge.cbs.dtu.dk/services/data.php
# Select resfinder.zip
# After extraction all data are in a folder named "resfinder"
# Alternaytive source
# https://bitbucket.org/genomicepidemiology/resfinder_db.git
# There are 15 classes of antibiotic resistance genes
#****************************************************************

#The database resfinder_db currently contains 15 classes of resistance genes
# Concatenate all fills into one single fasta file
# *** format resgenes database
# Skip these two lines if the data base has already been formatted.

# cat resfinder_db/*.fsa >resgenes.fasta
 makeblastdb -in resgenes.fasta -dbtype nucl

#*** blast search the combined database resgenes.fasta for the top hits

blastn -db resgenes.fasta -query $1 -out resgenes.blastout.txt -outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'

# Filter off identity ($3) <90% and coverage ($4/$10)<0.9
# This setting can be adjusted for the strigency of your search
awk '/Query/{qnam=$3;n=0}; !/#/&&$3>90&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10; print n,nam[n],res[n],id[n],rat[n]}' <resgenes.blastout.txt >tmp.1.txt

# Track of qstart and qend
# Filter off qstart and qend within 50 bps to the existing hits
# Theses settings can be adjusted for the strigency of your search
#awk '/Query/{qnam=$3;n=0}; !/#/&&$3>90&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10;flag=1; for(i=1;i<n;i++){delta1=(($8-qstart[i])^2)^0.5; delta2=(($9-qend[i])^2)^0.5; if(delta1<50||delta2<50){n--;flag=0;break}}if(flag==1){print n,nam[n],res[n],id[n],rat[n];flag=0}}' <resgenes.blastout.txt >tmp.1.txt


# Using flag to keep track of status, new query: flag=0, satisfying identities and coverage thrshold: flag=1, qstart, qend overlapping with existing hits: flag=2
awk '/Query/{qnam=$3;n=0;flag=0};  !/#/&&$3>90&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10; flag=1; for(i=1;i<n;i++){delta1=(($8-qstart[i])^2)^0.5; delta2=(($9-qend[i])^2)^0.5; if(delta1<50||delta2<50){n--;flag=2;break}}if(flag==1){print n,nam[n],res[n],id[n],rat[n];flag=2}}' <resgenes.blastout.txt >tmp.1.txt

# Consider 0 hits found,flag remains 0 after completing all hits
awk '/Query/{if(flag==0&&qnam!=""){print "0",qnam,"none"};qnam=$3;n=0;flag=0};  !/#/&&$3>90&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10; flag=1; for(i=1;i<n;i++){delta1=(($8-qstart[i])^2)^0.5; delta2=(($9-qend[i])^2)^0.5; if(delta1<50||delta2<50){n--;flag=2;break}}if(flag==1){print n,nam[n],res[n],id[n],rat[n];flag=2}}END{if(flag==0){print "0",qnam,"none"}}' <resgenes.blastout.txt >tmp.1.txt

# Consider 98% identities as a hit
# Compare tmp.1.txt at 90% vs tmp.1a.txt at 98%
awk '/Query/{if(flag==0&&qnam!=""){print "0",qnam,"none"};qnam=$3;n=0;flag=0};  !/#/&&$3>98&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10; flag=1; for(i=1;i<n;i++){delta1=(($8-qstart[i])^2)^0.5; delta2=(($9-qend[i])^2)^0.5; if(delta1<50||delta2<50){n--;flag=2;break}}if(flag==1){print n,nam[n],res[n],id[n],rat[n];flag=2}}END{if(flag==0){print "0",qnam,"none"}}' <resgenes.blastout.txt >tmp.1a.txt

# Consider 97% identities as a hit
# Compare tmp.1.txt at 90% vs tmp.1a.txt at 98% vs tmp.1b.txt at 97%
awk '/Query/{if(flag==0&&qnam!=""){print "0",qnam,"none"};qnam=$3;n=0;flag=0};  !/#/&&$3>97&&$4/$10>0.9{n++; nam[n]=$1; res[n]=$2; id[n]=$3; qstart[n]=$8; qend[n]=$9; rat[n]=$4/$10; flag=1; for(i=1;i<n;i++){delta1=(($8-qstart[i])^2)^0.5; delta2=(($9-qend[i])^2)^0.5; if(delta1<50||delta2<50){n--;flag=2;break}}if(flag==1){print n,nam[n],res[n],id[n],rat[n];flag=2}}END{if(flag==0){print "0",qnam,"none"}}' <resgenes.blastout.txt >tmp.1b.txt

# Listing all resistance genes in one single line
#awk '$2!=nam[n]{n++; nam[n]=$2; line[n]=nam[n]}; $2==nam[n]{line[n]=line[n]" "$3}END{for(i=1;i<=n;i++){print line[i]}}'<tmp.1.txt >tmp.2.txt

# Listing all resistance genes names in one single line
awk 'BEGIN{FS="[_]*[ ]*"}; $2!=nam[n]{n++; nam[n]=$2; line[n]=nam[n]}; $2==nam[n]{line[n]=line[n]" "$3}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.1.txt >tmp.2.txt

# Listing all resistance genes names in one single line at 98%
awk 'BEGIN{FS="[_]*[ ]*"}; $2!=nam[n]{n++; nam[n]=$2; line[n]=nam[n]}; $2==nam[n]{line[n]=line[n]" "$3}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.1a.txt >tmp.2a.txt

#DeLux version including counts
# Listing all resistance genes names in one single line use short names
awk 'BEGIN{FS="[_]*[ ]*"}; $2!=nam[n]{n++; nam[n]=$2; line[n]=""; ct[n]=0}; $2==nam[n]{ct[n]++;line[n]=line[n]" "$3}; END{for(i=1;i<=n;i++){print nam[i],ct[i]line[i]}}'<tmp.1.txt > final.txt

echo
cat final.txt
mv final.txt $1.result.txt
rm tmp*
