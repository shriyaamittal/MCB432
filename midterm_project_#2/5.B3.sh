head -n 100 all_otu_count.txt > all_otu_count_top_100.txt
awk '{print $2}' < all_otu_count_top_100.txt > otu_name_top_100.txt

awk -v nam="AF_A_8646_2_1" '$1==">"nam{print;exit}'< all.fasta 


