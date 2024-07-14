#!/usr/bin/env Rscript
library(vegan)

#load data
oe.tpm <- read.csv('../data/oe_tpm_counts_kallisto.csv', row.names = 1)
oe.tpm <- log(oe.tpm + 1)
dplex.tpm <- read.csv('../data/dpl_inf_tpm_counts_kallisto.csv', row.names = 1)
dplex.tpm <- log(dplex.tpm + 1)

#calculate distance matricies
#manhattan
oe.manhattan.matrix <- dist(oe.tpm, method='manhattan')
oe.manhattan.matrix <- as.matrix(oe.manhattan.matrix, labels=TRUE)
colnames(oe.manhattan.matrix) <- rownames(oe.manhattan.matrix) <- oe.tpm[['X']]

dplex.manhattan.matrix <- dist(dplex.tpm, method='manhattan')
dplex.manhattan.matrix <- as.matrix(dplex.manhattan.matrix, labels=TRUE)
colnames(dplex.manhattan.matrix) <- rownames(dplex.manhattan.matrix) <- dplex.tpm[['X']]

#reformat dataframes
oe.tpm <- t(oe.tpm)
dplex.tpm <- t(dplex.tpm)

#spearman
oe.spearman.matrix <- cor(oe.tpm, method="spearman")
oe.spearman.matrix <- (1 - oe.spearman.matrix)/2

dplex.spearman.matrix <- cor(dplex.tpm, method="spearman")
dplex.spearman.matrix <- (1 - dplex.spearman.matrix)/2

#pearson
oe.pearson.matrix <- cor(oe.tpm, method="pearson")
oe.pearson.matrix <- (1 - oe.pearson.matrix)/2

dplex.pearson.matrix <- cor(dplex.tpm, method="pearson")
dplex.pearson.matrix <- (1 - dplex.pearson.matrix)/2

#perform mantel tests
#manhattan
mantel.test.manhattan <- mantel(oe.manhattan.matrix, dplex.manhattan.matrix)
print(mantel.test.manhattan)

#spearman
mantel.test.spearman <- mantel(oe.spearman.matrix, dplex.spearman.matrix)
print(mantel.test.spearman)

#pearson
mantel.test.pearson <- mantel(oe.pearson.matrix, dplex.pearson.matrix)
print(mantel.test.pearson)

#Distance based-Redundancy analysis
# RDA
oe.pearson.dist <- as.dist(oe.pearson.matrix)
dplex.pearson.dist <- as.dist(dplex.pearson.matrix)
oe.pcoa <- cmdscale(oe.pearson.dist, eig = TRUE, k = nrow(oe.pearson.matrix) - 1)
dplex.pcoa <- cmdscale(dplex.pearson.dist, eig = TRUE, k = nrow(dplex.pearson.matrix) - 1)

# Extract the principal coordinates
oe.coords <- oe.pcoa$points
dplex.coords <- dplex.pcoa$points

#perform RDA
rda.model <- rda(oe.coords ~ dplex.coords)
summary(rda.model)

R2 <- RsquareAdj(rda.model)
print(R2)

#quick plots
oe.manhattan.vector <- as.vector(oe.manhattan.matrix)
dplex.manhattan.vector <- as.vector(dplex.manhattan.matrix)
plot(dplex.manhattan.vector,oe.manhattan.vector,
     main = "Mantel Test: Pairwise Distance Comparison",
     xlab = "OE",
     ylab = "Dplex",
     pch = 19,  # Plotting character
     col = "blue")  # Color of points

# Add a linear regression line
abline(lm(oe.manhattan.vector ~ dplex.manhattan.vector), col = "red")

oe.euclidian.vector <- as.vector(oe.euclidian.matrix)
dplex.euclidian.vector <- as.vector(dplex.euclidian.matrix)
plot(oe.euclidian.vector, dplex.euclidian.vector,
     main = "Mantel Test: Pairwise Distance Comparison",
     xlab = "OE",
     ylab = "Dplex",
     pch = 19,  # Plotting character
     col = "blue")  # Color of points

# Add a linear regression line
abline(lm(dplex.euclidian.vector ~ oe.euclidian.vector), col = "red")

oe.spearman.vector <- as.vector(oe.spearman.matrix)
dplex.spearman.vector <- as.vector(dplex.spearman.matrix)
spearman.df <- data.frame(oe.spearman.vector, dplex.spearman.vector)
write.csv(spearman.df, '../data/transcriptional_correlation_dplex_oe_spearman.csv', row.names = FALSE)
plot(oe.spearman.vector, dplex.spearman.vector,
     main = "Mantel Test: Pairwise Distance Comparison",
     xlab = "OE",
     ylab = "Dplex",
     pch = 19,  # Plotting character
     col = "blue")  # Color of points

# Add a linear regression line
abline(lm(dplex.spearman.vector ~ oe.spearman.vector), col = "red")

oe.pearson.vector <- as.vector(oe.pearson.matrix)
dplex.pearson.vector <- as.vector(dplex.pearson.matrix)
pearson.df <- data.frame(oe.pearson.vector, dplex.pearson.vector)
write.csv(pearson.df, '../data/transcriptional_correlation_dplex_oe_pearson.csv', row.names = FALSE)
plot(oe.pearson.vector, dplex.pearson.vector,
     main = "Mantel Test: Pairwise Distance Comparison",
     xlab = "OE",
     ylab = "Dplex",
     pch = 19,  # Plotting character
     col = "blue")  # Color of points

# Add a linear regression line
abline(lm(dplex.pearson.vector ~ oe.pearson.vector), col = "red")
