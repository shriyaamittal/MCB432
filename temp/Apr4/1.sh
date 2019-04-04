echo "Iteration 1"

## awk '/>/{n++}; (n<=100){print}; (n>100){exit}' <CA_A_all.fasta >seed_100.1.fasta 
awk '/>/{n++}; (n<=1000){print}; (n>1000){exit}' <CA_A_all.fasta >seed_1000.1.fasta 

## makeblastdb -in seed_100.1.fasta -dbtype nucl; blastn -query seed_100.1.fasta -db seed_100.1.fasta -out out_table.1.txt -outfmt '7 qseqid sseqid pident length qlen slen' 

makeblastdb -in seed_1000.1.fasta -dbtype nucl; blastn -query seed_1000.1.fasta -db seed_1000.1.fasta -out out_table.1.txt -outfmt '7 qseqid sseqid pident length qlen slen' 

awk '/Query/{n=0;print}; (!/#/&&$3>98&&$5/$4>=0.98){n++;print n, $0}' < out_table.1.txt > tmp.1.txt 

awk '/Query/{n++;nam[n]=$3}!/#/{ct[n]=$1; mem[n,$1]=$3}END{for(i=1;i<=n;i++){line=i" "nam[i]" "ct[i]; for(j=1;j<=ct[i];j++){line=line" "mem[i,j]}; print line}}'<tmp.1.txt >tmp.2.txt 

awk '{n=$1;nam[n]=$2;ct[n]=$3;line[n]=$0;for(i=1;i<=ct[n];i++){mem[n,i]=$(i+3)}for(j=1;j<n;j++){for(k=1;k<=ct[j];k++){if(mem[j,k]==nam[n]){if(ct[n]>=ct[j]){line[j]="";ct[j]=ct[n]};if(ct[n]<ct[j]){line[n]=""}}}}}END{for(l=1;l<=n;l++){if(line[l]!=""){print line[l]}}}'<tmp.2.txt >tmp.3.txt 

awk '{print $2,$3}'<tmp.3.txt| sort -k2nr|awk '{n++;print n, $1, $2, sum=sum+$2}'> cluster.1.txt 

rm new_seed.fasta; awk '{print $2}'<cluster.1.txt |while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' <seed_1000.1.fasta >>new_seed.1.fasta;done 

makeblastdb -in new_seed.1.fasta -dbtype nucl; blastn -query CA_A_all.fasta -db new_seed.1.fasta -out out_new.1.txt -max_target_seqs 1 -outfmt '7 qseqid sseqid pident length qlen slen' 

rm hit_list.1.txt nohit_list.1.txt; awk '!/#/{if($3>98&&$5/$4>=0.98){n++;print n, $0 >>"hit_list.1.txt"}else{m++; print m, $0 >>"nohit_list.1.txt"}}'< out_new.1.txt 

awk '{print $2}'<nohit_list.1.txt|while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' <CA_A_all.fasta >>rest.1.fasta;done 

echo "Iteration 2"

rm seed_1000.2.fasta;awk '{n++; $2}(n<=1000){print $2}(n>1000){exit}'<nohit_list.1.txt | while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' <rest.1.fasta >>seed_1000.2.fasta;done 

makeblastdb -in seed_1000.2.fasta -dbtype nucl; blastn -query seed_1000.2.fasta -db seed_1000.2.fasta -out out_table.2.txt -outfmt '7 qseqid sseqid pident length qlen slen' 

awk '/Query/{n=0;print}; (!/#/&&$3>98&&$5/$4>=0.98){n++;print n, $0}' <out_table.2.txt| awk '/Query/{n++;nam[n]=$3}!/#/{ct[n]=$1; mem[n,$1]=$3}END{for(i=1;i<=n;i++){line[i]=i" "nam[i]" "ct[i]; for(j=1;j<=ct[i];j++){line[i]=line[i]" "mem[i,j]};for(j=1;j<i;j++){for(k=1;k<=ct[j];k++){if(mem[j,k]==nam[i]){if(ct[i]>=ct[j]){line[j]="";ct[j]=ct[i]};if(ct[i]<ct[j]){line[i]=""}}}}}for(i=1;i<=n;i++){if(line[i]!=""){print line[i]}}}'| awk '{print $2,$3}'| sort -k2nr|awk '{n++;print n, $1, $2, sum=sum+$2}'>cluster.2.txt 

rm new_seed.2.fasta; awk '{print $2}'<cluster.2.txt |while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' <seed_1000.2.fasta >>new_seed.2.fasta;done 

makeblastdb -in new_seed.2.fasta -dbtype nucl; blastn -query rest.1.fasta -db new_seed.2.fasta -out out_new.2.txt -max_target_seqs 1 -outfmt '7 qseqid sseqid pident length qlen slen' 

echo "Iteration 3"

awk '{print $2}'<hit_list.1.txt|while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' <CA_A_all.fasta >>rest.2.fasta;done

cp rest.2.fasta seed_1000.3.fasta; makeblastdb -in seed_1000.3.fasta -dbtype nucl; blastn -query seed_1000.3.fasta -db seed_1000.3.fasta -out out_table.3.txt -outfmt '7 qseqid sseqid pident length qlen slen' 

awk '/Query/{n=0;print}; (!/#/&&$3>98&&$5/$4>=0.98){n++;print n, $0}' <out_table.3.txt| awk '/Query/{n++;nam[n]=$3}!/#/{ct[n]=$1; mem[n,$1]=$3}END{for(i=1;i<=n;i++){line[i]=i" "nam[i]" "ct[i]; for(j=1;j<=ct[i];j++){line[i]=line[i]" "mem[i,j]};for(j=1;j<i;j++){for(k=1;k<=ct[j];k++){if(mem[j,k]==nam[i]){if(ct[i]>=ct[j]){line[j]="";ct[j]=ct[i]};if(ct[i]<ct[j]){line[i]=""}}}}}for(i=1;i<=n;i++){if(line[i]!=""){print line[i]}}}'| awk '{print $2,$3}'| sort -k2nr|awk '{n++;print n, $1, $2, sum=sum+$2}'>cluster.3.txt 

rm new_seed.3.fasta; awk '{print $2}'<cluster.3.txt |while read p; do awk -v nam=$p 'BEGIN{RS=">"}($1==nam){print ">"$0;exit}' < seed_1000.3.fasta >>new_seed.3.fasta;done 

cat new_seed.1.fasta new_seed.2.fasta new_seed.3.fasta>>all_seed.fasta 

rm all.cluster.txt;cat cluster.1.txt cluster.2.txt cluster.3.txt>>all.cluster.txt; sort -k3nr all.cluster.txt|awk '{n++; print n, $2,$3,sum=sum+$3}'>final.tx


