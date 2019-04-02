library(vegan)
library(MASS)

data(varespec)
vare.dis <- vegdist(varespec)
vare.mds0 <- isoMDS(vare.dis)

png("figure1.png")
stressplot(vare.mds0,vare.dis)
dev.off()


png("figure2.png")
ordiplot(vare.mds0, type = "t")
dev.off()

vare.mds <- metaMDS(varespec, trace = FALSE)

vare.mds

png("figure3.png")
plot(vare.mds, type = "t")
dev.off()

vare.pca <- rda(varespec)
vare.pca

png('figure4.png')
plot(vare.pca)
dev.off()

sum(apply(varespec, 2, var))

png('figure5.png')
biplot(vare.pca, scaling = -1)
dev.off()

vare.pca <- rda(varespec, scale = TRUE)
vare.pca

png('figure6.png')
plot(vare.pca, scaling = 3)
dev.off()

dim(varespec)

vare.ca <- cca(varespec)
vare.ca

png('figure7.png')
plot(vare.ca)
dev.off()

chisq.test(varespec/sum(varespec))

png('figure8.png')
plot(vare.ca, scaling = 1)
dev.off()

data(varechem)
ef <- envfit(vare.mds, varechem, permu = 999)
ef

png('figure9.png')
plot(vare.mds, display = "sites")
dev.off()

data(dune)
data(dune.env)
dune.ca <- cca(dune)

ef <- envfit(dune.ca, dune.env, permutations = 999)
ef

png('figure10.png')
plot(dune.ca, display = "sites")
dev.off()


png('figure11.png')
plot(dune.ca, display = "sites", type = "p")
with(dune.env, ordiellipse(dune.ca, Management, kind = "se", conf = 0.95))
with(dune.env, ordispider(dune.ca, Management, col = "blue", label= TRUE))
with(dune.env, ordihull(dune.ca, Management, col="blue", lty=2))
dev.off()

