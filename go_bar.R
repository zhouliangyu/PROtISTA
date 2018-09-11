# temporary code
# ==============

df <- read.table("./scd49_ce_go.txt", sep="\t", skip=10, header=TRUE)
df <- df[c(1,6)]
df <- df[1:20,]
colnames(df) <- c("GO_TERM", "Fold_enrichment")
tempMar <- par()$mar
dev.off()
pdf("go_bar.pdf", height=5, width=7)
par(mar = c(5, 18, 4, 3))
barplot(rev(as.numeric(as.character(df$Fold_enrichment))), names.arg=rev(df$GO_TERM), horiz=TRUE,
	main="C. elegans, S/TQ interval: 49 AA", xlim=c(0, 12), las=1,
	cex.names=0.8, xlab="Fold Enrichment")
par(mar =  tempMar)
dev.off()
