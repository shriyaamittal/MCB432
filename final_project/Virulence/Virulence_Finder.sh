makeblastdb -in VFDB_setA_nt.fas -dbtype nucl

blastn -db VFDB_setA_nt.fas -query $1 -out virulence.blastout.txt -num_descriptions 5 -num_alignments 1
##-outfmt '7 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore'

awk -v nam=$p '/Query=/{Qu=$2}/>/{print Qu,$2}/No hits found/{print Qu,nam"_0"}' <virulence.blastout.txt>virulence.hit.txt

cp virulence.hit.txt tmp.1.txt

awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.1.txt >tmp.2.txt

awk 'BEGIN{FS="[ ]*[_]*"}{print $1,$3,$5,$7,$9,$11,$13,$15}' <tmp.2.txt >tmp.3.txt

awk '$1!=nam[n]{n++; nam[n]=$1;line[n]=$1}; $1==nam[n]{line[n]=line[n]" "$2}; END{for(i=1;i<=n;i++){print line[i]}}'<tmp.1.txt| awk 'BEGIN{FS="[ ]*[_]*"}{print $1"\t"$3"\t"$5"\t"$7"\t"$9"\t"$11"\t"$13"\t"$15}'>tmp.4.txt

cat tmp.4.txt | tr -d ')' > final.txt

echo 
cat final.txt
mv final.txt $1.result.txt
rm tmp*
