cat *.fasta > all.fasta

awk 'BEGIN{RS=">"}$2==""{print}'< all.fasta
