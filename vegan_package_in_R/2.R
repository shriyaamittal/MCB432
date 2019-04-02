## After finishing midterm project #2 analysis
## use the following

plot(OTU.rda, display = "sites", type = "p")
with(x2, ordiellipse(OTU.rda, Color, kind = "se", conf = 0.95))
with(x2, ordispider(OTU.rda, Color, col = "blue", label= TRUE))
with(x2, ordihull(OTU.rda, Color, col="blue", lty=2))


plot(OTU.rda, display = "sites", type = "p")
with(x2, ordiellipse(OTU.rda, Symbol, kind = "se", conf = 0.95))
with(x2, ordispider(OTU.rda, Symbol, col = "blue", label= TRUE))
with(x2, ordihull(OTU.rda, Symbol, col="blue", lty=2))
