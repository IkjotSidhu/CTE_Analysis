library(tidyverse)
library(RobustRankAggreg)

PD_GSEA <- read_csv("PD_GSEA.csv")
HD_GSEA <- read_csv("HD_GSEA.csv")
HD_DE <- read_csv("HD_DE.csv") %>% filter(!is.na(symbol))
#HD_GSEA_list <- as.vector(HD_GSEA$leadingEdge)
HD_glist <- str_split(HD_GSEA$leadingEdge,pattern=" ")
PD_glist <- str_split(PD_GSEA$leadingEdge,pattern=" ")
CTES_RHIN_glist <- str_split(CTES_RHIN_GSEA$leadingEdge,pattern=" ")

#PD_DE <- read_csv("PD_DE.csv") %>% filter(!is.na(symbol))
CTES_RHIN_GSEA <- read_csv("CTES_RHIN_all_filthi_de_annot_gsea_c2cp.csv")

ranked_HD <- aggregateRanks(HD_glist)
ranked_PD <- aggregateRanks(PD_glist)
ranked_CTES_GSEA <- aggregateRanks(CTES_RHIN_glist)


# Make sample input data
glist <- list(sample(letters, 4), sample(letters, 10), sample(letters, 12))
# Aggregate the inputs
aggregateRanks(glist = glist, N = length(letters))
aggregateRanks(glist = glist, N = length(letters), method = "stuart")
# Since we know the cutoffs for the lists in advance (4, 10, 12) we can use
# the more accurate algorithm with parameter topCutoff
# Use the rank matrix instead of the gene lists as the input
r = rankMatrix(glist)