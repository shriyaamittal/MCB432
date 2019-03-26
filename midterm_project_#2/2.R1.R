library(ape)
library(kmer)

all.dna<- read.dna("all.fasta",format="fasta")

#Sys.time();all.otu<- otu(all.dna, k=5, threshold=0.97, method="farthest", nstart=20);Sys.time()
Sys.time();all.cl<- cluster(all.dna, k=5, nstart=20);Sys.time()
Sys.time();set.seed(999);all.otu<- otu(all.dna, k=5, threshold=0.97, method="farthest", nstart=20);Sys.time()
Sys.time();set.seed(999);all.tree<- cluster(all.dna, k=5, nstart=20);Sys.time()

png('figure1.png')
Sys.time();plot(all.tree,main="all.tree by kmer cluster",horiz=T);Sys.time()
dev.off()

all.tre<- as.phylo(all.cl)

write.nexus(all.tre, file="all.nex")
write.table(all.otu,"all_otu.txt")


