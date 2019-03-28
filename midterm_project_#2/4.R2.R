all_otu_count<- read.table("all_otu_count.txt")
colnames(all_otu_count)<-c("OTU","name","all","ct")

otu.1.count<- read.table("otu.AF_A_2.count.txt")
colnames(otu.1.count)<-c("OTU","AF_A","ct")

otu.2.count<- read.table("otu.AF_B_2.count.txt")
colnames(otu.2.count)<-c("OTU","AF_B","ct")

otu.3.count<- read.table("otu.AM_A_4.count.txt")
colnames(otu.3.count)<-c("OTU","AM_A","ct")

otu.4.count<- read.table("otu.AM_B_4.count.txt")
colnames(otu.4.count)<-c("OTU","AM_B","ct")

otu.5.count<- read.table("otu.AQ_A_1.count.txt")
colnames(otu.5.count)<-c("OTU","AQ_A","ct")

otu.6.count<- read.table("otu.AQ_B_1.count.txt")
colnames(otu.6.count)<-c("OTU","AQ_B","ct")

otu.7.count<- read.table("otu.AX_A_2.count.txt")
colnames(otu.7.count)<-c("OTU","AX_A","ct")

otu.8.count<- read.table("otu.AX_B_2.count.txt")
colnames(otu.8.count)<-c("OTU","AX_B","ct")

otu.9.count<- read.table("otu.BB_A_3.count.txt")
colnames(otu.9.count)<-c("OTU","BB_A","ct")

otu.10.count<- read.table("otu.BB_B_3.count.txt")
colnames(otu.10.count)<-c("OTU","BB_B","ct")

otu.11.count<- read.table("otu.BF_A_3.count.txt")
colnames(otu.11.count)<-c("OTU","BF_A","ct")

otu.12.count<- read.table("otu.BF_B_3.count.txt")
colnames(otu.12.count)<-c("OTU","BF_B","ct")

otu.13.count<- read.table("otu.BI_A_3.count.txt")
colnames(otu.13.count)<-c("OTU","BI_A","ct")

otu.14.count<- read.table("otu.BI_B_3.count.txt")
colnames(otu.14.count)<-c("OTU","BI_B","ct")

otu.15.count<- read.table("otu.BM_A_4.count.txt")
colnames(otu.15.count)<-c("OTU","BM_A","ct")

otu.16.count<- read.table("otu.BM_B_4.count.txt")
colnames(otu.16.count)<-c("OTU","BM_B","ct")

otu.17.count<- read.table("otu.BU_A_4.count.txt")
colnames(otu.17.count)<-c("OTU","BU_A","ct")

otu.18.count<- read.table("otu.BU_B_4.count.txt")
colnames(otu.18.count)<-c("OTU","BU_B","ct")

otu.19.count<- read.table("otu.BX_A_3.count.txt")
colnames(otu.19.count)<-c("OTU","BX_A","ct")

otu.20.count<- read.table("otu.BX_B_3.count.txt")
colnames(otu.20.count)<-c("OTU","BX_B","ct")

m1<-merge(x=all_otu_count,y=otu.1.count,by="OTU", all.x=T, all.y=T)
m1[is.na(m1)]<-0; m1[,6]<-NULL; m1[,4]<-NULL; head(m1)

otu.2.count[,3]<-NULL
otu.3.count[,3]<-NULL
otu.4.count[,3]<-NULL
otu.5.count[,3]<-NULL
otu.6.count[,3]<-NULL
otu.7.count[,3]<-NULL
otu.8.count[,3]<-NULL
otu.9.count[,3]<-NULL
otu.10.count[,3]<-NULL
otu.11.count[,3]<-NULL
otu.12.count[,3]<-NULL
otu.13.count[,3]<-NULL
otu.14.count[,3]<-NULL
otu.15.count[,3]<-NULL
otu.16.count[,3]<-NULL
otu.17.count[,3]<-NULL
otu.18.count[,3]<-NULL
otu.19.count[,3]<-NULL
otu.20.count[,3]<-NULL

m2<-merge(x=m1,y=otu.2.count,by="OTU", all.x=T, all.y=T);m2[is.na(m2)]<-0;head(m2)
m3<-merge(x=m2,y=otu.3.count,by="OTU", all.x=T, all.y=T);m3[is.na(m3)]<-0;head(m3)
m4<-merge(x=m3,y=otu.4.count,by="OTU", all.x=T, all.y=T);m4[is.na(m4)]<-0;head(m4)
m5<-merge(x=m4,y=otu.5.count,by="OTU", all.x=T, all.y=T);m5[is.na(m5)]<-0;head(m5)
m6<-merge(x=m5,y=otu.6.count,by="OTU", all.x=T, all.y=T);m6[is.na(m6)]<-0;head(m6)
m7<-merge(x=m6,y=otu.7.count,by="OTU", all.x=T, all.y=T);m7[is.na(m7)]<-0;head(m7)
m8<-merge(x=m7,y=otu.8.count,by="OTU", all.x=T, all.y=T);m8[is.na(m8)]<-0;head(m8)
m9<-merge(x=m8,y=otu.9.count,by="OTU", all.x=T, all.y=T);m9[is.na(m9)]<-0;head(m9)
m10<-merge(x=m9,y=otu.10.count,by="OTU", all.x=T, all.y=T);m10[is.na(m10)]<-0;head(m10)
m11<-merge(x=m10,y=otu.11.count,by="OTU", all.x=T, all.y=T);m11[is.na(m11)]<-0;head(m11)
m12<-merge(x=m11,y=otu.12.count,by="OTU", all.x=T, all.y=T);m12[is.na(m12)]<-0;head(m12)
m13<-merge(x=m12,y=otu.13.count,by="OTU", all.x=T, all.y=T);m13[is.na(m13)]<-0;head(m13)
m14<-merge(x=m13,y=otu.14.count,by="OTU", all.x=T, all.y=T);m14[is.na(m14)]<-0;head(m14)
m15<-merge(x=m14,y=otu.15.count,by="OTU", all.x=T, all.y=T);m15[is.na(m15)]<-0;head(m15)
m16<-merge(x=m15,y=otu.16.count,by="OTU", all.x=T, all.y=T);m16[is.na(m16)]<-0;head(m16)
m17<-merge(x=m16,y=otu.17.count,by="OTU", all.x=T, all.y=T);m17[is.na(m17)]<-0;head(m17)
m18<-merge(x=m17,y=otu.18.count,by="OTU", all.x=T, all.y=T);m18[is.na(m18)]<-0;head(m18)
m19<-merge(x=m18,y=otu.19.count,by="OTU", all.x=T, all.y=T);m19[is.na(m19)]<-0;head(m19)
m20<-merge(x=m19,y=otu.20.count,by="OTU", all.x=T, all.y=T);m20[is.na(m20)]<-0;head(m20)

write.table(m20,"otu_table.txt")

OTU<-m20[,-2]
OTU_table<-OTU[,-1]
OTU_tablex<-OTU_table[,-1]; head(OTU_tablex)
write.csv(OTU_tablex,"otu_tablex.csv")

library(ape)
library(vegan)

OTUx<-t(OTU_tablex) # Transpose

##PCoA
OTU.dist<- vegdist(OTUx)
Sys.time();OTU.pcoa<- pcoa(OTU.dist);Sys.time();
OTU.pcoa

png('figure2.png')#,res=500,width = 3, height = 3, units = 'in')
biplot(OTU.pcoa)
#points=c(16, 1)
dev.off()

## RDA
OTU.rda<- rda(OTUx,scale=TRUE)

x2<- read.csv("meta01.csv",row.names=1)

envfit(OTU.rda,x2[,1:2],perm = 999)

png('figure3.png')
biplot(OTU.rda,scaling=-3,type=c("t","p"),display="spe",col="gray", main="PCA analysis for 8 Sets",sub="keys= X(Red) Y(Green) Z(Blue) o(Circle) p(Diamond) q(Triangle)")
points(OTU.rda,scaling=-3,pch=x2$Symbol,bg=x2$Color,cex=1.5)
text(OTU.rda,scaling=-3,pos=1,cex=1.2)
dev.off()

x2<- read.csv("meta02.csv",row.names=1)

envfit(OTU.rda,x2[,1:2],perm = 999)

png('figure4.png')
biplot(OTU.rda,scaling=-3,type=c("t","p"),display="spe",col="gray", main="PCA analysis for 8 Sets",sub="keys= A(Black) B(Red) o(Circle) p(Diamond) q(Triangle)")
points(OTU.rda,scaling=-3,pch=x2$Symbol,bg=x2$Color,cex=1.5)
text(OTU.rda,scaling=-3,pos=1,cex=1.2)
dev.off()


x2<- read.csv("meta03.csv",row.names=1)

envfit(OTU.rda,x2[,1:2],perm = 999)

png('figure5.png')
biplot(OTU.rda,scaling=-3,type=c("t","p"),display="spe",col="gray", main="PCA analysis for 8 Sets",sub="keys= X(Red) Y(Green) Z(Blue) A(Circle) B(Triangle)")
points(OTU.rda,scaling=-3,pch=x2$Symbol,bg=x2$Color,cex=1.5)
text(OTU.rda,scaling=-3,pos=1,cex=1.2)
dev.off()

