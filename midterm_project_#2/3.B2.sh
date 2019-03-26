awk 'BEGIN{FS="[\:]*[ ]*";gsub(/\"/,"")}$(NF)>max{max=$(NF)}/\*/{nam[$(NF)]=substr($1,2,length($1)-2)}{ct[$(NF)]++} END{for(i=1;i<=max;i++){print i,nam[i],ct[i],sum=sum+ct[i]}}' < all_otu.txt >all_otu_count.txt

awk 'BEGIN{FS="[\:]*[ ]*"}$(NF)>max{max=$(NF)}/\*/{nam[$(NF)]=substr($1,2,length($1)-2)}{ct[$(NF)]++} END{for(i=1;i<=max;i++){print i,nam[i],ct[i],sum=sum+ct[i]}}' < all_otu.txt|sort -k3nr >all_otu_count.txt
## This file is used is to plot OTU abundance.



