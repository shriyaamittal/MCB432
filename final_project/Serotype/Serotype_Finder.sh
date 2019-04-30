#!/bin/sh

#****************************************************************
#  serotype_finder.2.sh
#
#  Created by Meng Ho on 4/9/19.
#  This script is intended for the MCB432 Course as an
#  In-class shell-scripting exercise
#  Contact mho1@illinois.edu for comments.
#
# Input file is a fasta file containing multiple sequences
# Usage: serotype_finder.2.sh file.fasta
# Download database from https://cge.cbs.dtu.dk/services/data.php
# Select serotypefinder.zip
#****************************************************************

#*** You may need to modify the script for your database name
#*** format serotype database
#Skip these two lines if the data base has already been formatted.
#makeblastdb -in O_antigens.fsa -dbtype nucl
#makeblastdb -in H_antigens.fsa -dbtype nucl
makeblastdb -in O_type.fsa -dbtype nucl
makeblastdb -in H_type.fsa -dbtype nucl

#*** blast search each Antigen dataset for the top hits
#rm H-antigens.blastout.txt O-antigens.blastout.txt
#blastn -db H-antigens.fsa -query $1 -out H-antigens.blastout.txt -num_descriptions 5 -num_alignments 5
#blastn -db O-antigens.fsa -query $1 -out O-antigens.blastout.txt -num_descriptions 20 -num_alignments 20
blastn -db H_type.fsa -query $1 -out H-antigens.blastout.txt -num_descriptions 5 -num_alignments 5
blastn -db O_type.fsa -query $1 -out O-antigens.blastout.txt -num_descriptions 20 -num_alignments 20
awk '/Query= /{print}/>/{print $2,substr($2,4,length($2)-3)}/Length/{print}/Identities =/{print}'<H-antigens.blastout.txt >H-serotype_result.txt
awk '/Query= /{print}/>/{print $2,substr($2,4,length($2)-3)}/Length/{print}/Identities =/{print}'<O-antigens.blastout.txt >O-serotype_result.txt


#rm H-antigens.blastout.txt O-antigens.blastout.txt

# Run blastn and output tabular for H-antigen
#blastn -db H-antigens.fsa -query $1 -out H-antigens.blastout.txt -max_target_seqs 1 -outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'
blastn -db H_type.fsa -query $1 -out H-antigens.blastout.txt -max_target_seqs 1 -outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'
#Parse H-antigen blastn output table
awk '/Query/{nam=$3}!/#/{gsub(/_/," ");print $1,$2,$5,$6,$7/$13}'<H-antigens.blastout.txt >tmp.1.txt

# Run blastn and output tabular for O-antigen
#blastn -db O-antigens.fsa -query $1 -out O-antigens.blastout.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'
blastn -db O_type.fsa -query $1 -out O-antigens.blastout.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'

# Parse O-antigen blastn output table and merge with H-antigen result
awk '/Query/{nam=$3;n=0}!/#/{gsub(/_/," ");if($2!=hit[n]&&$2!=hit[n-1]){n++;print $1,hit[n]=$2,$5,$6,$7/$13}}'<O-antigens.blastout.txt >>tmp.1.txt
sort tmp.1.txt >tmp.2.txt

# Writing antigen determinant genes into single line
awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]"\t"$2"_"$3}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.2.txt >tmp.3.txt

# Writing a formatted output
awk 'BEGIN{printf "%-15s %-15s %-15s\n", "Accession","H-Antigen_type","O-Antigen_type"}{printf "%-15s %-15s %-15s %-15s\n", $1,$2,$3,$4}'  <tmp.3.txt > final.txt

echo
cat final.txt
mv final.txt $1.result.txt
rm tmp*

