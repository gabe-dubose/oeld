library(dplyr)

# function to extract entries based on sample id
select.stage <- function(matrix, char) {
  entries <- rownames(matrix)
  subset(entries, substr(entries, 6, 6) == char)
}

# function to get average distance
distance.between.stages <- function(matrix, stage1, stage2){
  stage1.entries <- select.stage(matrix, stage1)
  stage2.entries <- select.stage(matrix, stage2)
  average.distance <- mean(matrix[stage1.entries, stage2.entries])
  return(average.distance)
}

#load data
oe.tpm <- read.csv('../data/oe_tpm_counts_kallisto.csv', row.names = 1)
dplex.tpm <- read.csv('../data/dpl_inf_tpm_counts_kallisto.csv', row.names = 1)

#calculate manhattan distances
oe.manhattan.matrix <- dist(oe.tpm, method='manhattan')
oe.manhattan.matrix <- as.matrix(oe.manhattan.matrix, labels=TRUE)

dplex.manhattan.matrix <- dist(dplex.tpm, method='manhattan')
dplex.manhattan.matrix <- as.matrix(dplex.manhattan.matrix, labels=TRUE)

#get distances OE
distance.between.stages(oe.manhattan.matrix, '3', '5')
distance.between.stages(oe.manhattan.matrix, '5', 'E')
distance.between.stages(oe.manhattan.matrix, 'E', 'L')
distance.between.stages(oe.manhattan.matrix, 'L', 'A')

#get dplex distances
distance.between.stages(dplex.manhattan.matrix, '3', '5')
distance.between.stages(dplex.manhattan.matrix, '5', 'E')
distance.between.stages(dplex.manhattan.matrix, 'E', 'L')
distance.between.stages(dplex.manhattan.matrix, 'L', 'A')

#log transform
oe.tpm <- log(oe.tpm + 1)
dplex.tpm <- log(dplex.tpm + 1)

#reformat dataframes
oe.tpm <- t(oe.tpm)
dplex.tpm <- t(dplex.tpm)

#pearson
oe.pearson.matrix <- cor(oe.tpm, method="pearson")
oe.pearson.matrix <- (1 - oe.pearson.matrix)/2

dplex.pearson.matrix <- cor(dplex.tpm, method="pearson")
dplex.pearson.matrix <- (1 - dplex.pearson.matrix)/2

#get distances OE
distance.between.stages(oe.pearson.matrix, '3', '5')
distance.between.stages(oe.pearson.matrix, '5', 'E')
distance.between.stages(oe.pearson.matrix, 'E', 'L')
distance.between.stages(oe.pearson.matrix, 'L', 'A')

#get dplex distances
distance.between.stages(dplex.pearson.matrix, '3', '5')
distance.between.stages(dplex.pearson.matrix, '5', 'E')
distance.between.stages(dplex.pearson.matrix, 'E', 'L')
distance.between.stages(dplex.pearson.matrix, 'L', 'A')

#spearman
oe.spearman.matrix <- cor(oe.tpm, method="spearman")
oe.spearman.matrix <- (1 - oe.spearman.matrix)/2

dplex.spearman.matrix <- cor(dplex.tpm, method="spearman")
dplex.spearman.matrix <- (1 - dplex.spearman.matrix)/2

#get distances OE
distance.between.stages(oe.spearman.matrix, '3', '5')
distance.between.stages(oe.spearman.matrix, '5', 'E')
distance.between.stages(oe.spearman.matrix, 'E', 'L')
distance.between.stages(oe.spearman.matrix, 'L', 'A')

#get dplex distances
distance.between.stages(dplex.spearman.matrix, '3', '5')
distance.between.stages(dplex.spearman.matrix, '5', 'E')
distance.between.stages(dplex.spearman.matrix, 'E', 'L')
distance.between.stages(dplex.spearman.matrix, 'L', 'A')
