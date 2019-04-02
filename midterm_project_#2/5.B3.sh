head -n 100 all_otu_count.txt > all_otu_count_top_100.txt
awk '{print $2}' < all_otu_count_top_100.txt > otu_name_top_100.txt


#awk -v nam="AF_A_8646_2_1" '$1==">"nam{print;exit}'< all.fasta 
#awk -v nam="AF_A_8646_2_1" 'BEGIN{RS=">"} $1==nam{print $0;exit}'<all.fasta 
#awk -v nam="AF_A_8646_2_1*" 'BEGIN{RS=">"} $1==substr(nam,1,length(nam)-1){print $0;exit}'<all.fasta 

# The out.fasta file must have 100 sequences
rm out.fasta
cat otu_name_top_100.txt | while read p; do awk -v nam=$p 'BEGIN{RS=">"} $1==substr(nam,1,length(nam)-1){print ">"$0;exit}'<all.fasta >> out.fasta ; done

blastn -p 
makeblastdb -in out.fasta -dbtype nucl
blastn -query out.fasta -db out.fasta -outfmt 7 -out blasttable.txt


