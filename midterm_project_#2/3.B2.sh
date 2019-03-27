#awk 'BEGIN{FS="[\:]*[ ]*";gsub(/\"/,"")}$(NF)>max{max=$(NF)}/\*/{nam[$(NF)]=substr($1,2,length($1)-2)}{ct[$(NF)]++} END{for(i=1;i<=max;i++){print i,nam[i],ct[i],sum=sum+ct[i]}}' < all_otu.txt >all_otu_count.txt

#awk 'BEGIN{FS="[\:]*[ ]*"}$(NF)>max{max=$(NF)}/\*/{nam[$(NF)]=substr($1,2,length($1)-2)}{ct[$(NF)]++} END{for(i=1;i<=max;i++){print i,nam[i],ct[i],sum=sum+ct[i]}}' < all_otu.txt|sort -k3nr >all_otu_count.txt
## This file is used is to plot OTU abundance.

## Create split txt files
#while read line;
#do
#	arrIN=(${line//./ })
#	echo $arrIN
#	echo "otu."$arrIN".txt"
#	touch "otu."$arrIN".txt"
#done<listFasta

## Populate split txt files
#while read line;
#do
#	arrIN=(${line// / })
#	arrIN2=(${arrIN//\"/ })
#	arrIN3=(${arrIN2//_/ })
#	fnam="otu."${arrIN3[0]}"_"${arrIN3[1]}"_"${arrIN3[3]}".txt"
#	echo $line >> ${fnam}
#done<all_otu.txt

ls otu.*.txt|awk '{gsub(/\.txt/,"");print}'|while read p; do echo $p;awk 'BEGIN{FS="[\:]*[ ]*"}$(NF)>max{max=$(NF)}/\*/{nam[$(NF)]=substr($1,2,length($1)-2)}{ct[$(NF)]++}END{for(i=1;i<=max;i++){if(ct[i]>0){print i,ct[i],sum=sum+ct[i]}}}' < $p.txt|sort -k2nr >$p.count.txt;done

rm otu.x__.count.txt otu.x__.txt x__.txt


