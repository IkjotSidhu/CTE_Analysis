library(tidyverse)
library(RobustRankAggreg)
library(pheatmap)

PD_GSEA <- read_csv("PD_GSEA.csv")
HD_GSEA <- read_csv("HD_GSEA.csv")

# Import DE results from HD,PD and CTES vs RHIN
HD_DE <- read_csv("HD_DE.csv") %>% filter(!is.na(symbol))
PD_DE <- read_csv("PD_DE.csv") %>% filter(!is.na(symbol))
CTES_RHIN_DE <- read_csv("CTES_RHIN_all_filthi_de_annot.csv") %>% filter(!is.na(gene_name))
CTES_RHIN_DE <-  arrange(CTES_RHIN_DE,Status__CTES__pvalue) %>% mutate(CTE_RHIN_Rank=row_number())

CTEM_RHIN_DE <- read_csv("CTEM_RHIN_all_filthi_de_annot.csv") %>% filter(!is.na(gene_name))
CTEM_RHIN_DE <-  arrange(CTEM_RHIN_DE,Status__CTEM__pvalue) %>% mutate(CTE_RHIN_Rank=row_number())

CTEM_RHIN_GSEA <- read_csv("CTEM_RHIN_all_filthi_de_annot_gsea_c2cp.csv") %>% arrange(pval)


#HD_GSEA_list <- as.vector(HD_GSEA$leadingEdge)

# Import GSEA pathways
HD_glist <- str_split(HD_GSEA$leadingEdge,pattern=" ")
PD_glist <- str_split(PD_GSEA$leadingEdge,pattern=" ")
CTES_RHIN_glist <- str_split(CTES_RHIN_GSEA$leadingEdge,pattern=" ")

#PD_DE <- read_csv("PD_DE.csv") %>% filter(!is.na(symbol))
CTES_RHIN_GSEA <- read_csv("CTES_RHIN_all_filthi_de_annot_gsea_c2cp.csv") %>% arrange(pval)

ranked_HD <- aggregateRanks(HD_glist)
ranked_PD <- aggregateRanks(PD_glist)
ranked_CTES_GSEA <- aggregateRanks(CTES_RHIN_glist)

HD_pathways <- (HD_GSEA$pathway)
PD_pathways <- (PD_GSEA$pathway)
CTE_RHIN_pathways <- (CTES_RHIN_GSEA$pathway)


glist <- append(HD_glist,PD_glist)
g1_list <- append(glist,CTES_RHIN_glist)
g2_list <- list(HD_pathways,PD_pathways,CTE_RHIN_pathways)

ranked_pathway_list <- aggregateRanks(g2_list)

# Running RRAGRR
HD_pathways_sorted <- (fgseaResTidy_HD$pathway)
PD_pathways_sorted <- (fgseaResTidy_PD$pathway)
CTES_pathway_sorted <- (CTES_RHIN_GSEA$pathway)
CTEM_pathway_sorted <- (CTEM_RHIN_GSEA$pathway)

g3_list <- list(HD_pathways_sorted,PD_pathways_sorted,CTEM_pathway_sorted,CTES_pathway_sorted)

ranked_pathway_list_4 <- aggregateRanks(g3_list)

ranked_pathway_list_4 <- remove_rownames(ranked_pathway_list_4)

RRA_result_GSEA <- ranked_pathway_list_4 %>% arrange(Score)

# ADDING COLUMNS
# Matching NES from HD
RRA_result_GSEA$NES_HD <- fgseaResTidy_HD$NES[match(RRA_result_GSEA$Name, fgseaResTidy_HD$pathway)]

# Matching NES from PD
RRA_result_GSEA$NES_PD <- fgseaResTidy_PD$NES[match(RRA_result_GSEA$Name, fgseaResTidy_PD$pathway)]

# Matching NES from CTES
RRA_result_GSEA$NES_CTES <- CTES_RHIN_GSEA$NES[match(RRA_result_GSEA$Name, CTES_RHIN_GSEA$pathway)]

# Matching NES from CTEM
RRA_result_GSEA$NES_CTEM <- CTEM_RHIN_GSEA$NES[match(RRA_result_GSEA$Name, CTEM_RHIN_GSEA$pathway)]

# Matching pval from HD
RRA_result_GSEA$padj_HD <- fgseaResTidy_HD$padj[match(RRA_result_GSEA$Name, fgseaResTidy_HD$pathway)]

# Matching pval from PD
RRA_result_GSEA$padj_PD <- fgseaResTidy_PD$padj[match(RRA_result_GSEA$Name, fgseaResTidy_PD$pathway)]

# Matching pval from CTEM
RRA_result_GSEA$padj_CTEM <- CTEM_RHIN_GSEA$padj[match(RRA_result_GSEA$Name, CTEM_RHIN_GSEA$pathway)]

# Matching pval from CTES
RRA_result_GSEA$padj_CTES <- CTES_RHIN_GSEA$padj[match(RRA_result_GSEA$Name, CTES_RHIN_GSEA$pathway)]

write_csv(RRA_result_GSEA,"RRA_result_mapped.csv")

# Cleaning up the data
df %>% mutate(g = case_when(a == 2 | a == 5 | a == 7 | (a == 1 & b == 4) ~ 2,
                            a == 0 | a == 1 | a == 4 | a == 3 |  c == 4 ~ 3,
                            TRUE ~ NA_real_))

RRA_result_GSEA %>% mutate(NES_HD=case_when(padj_HD > 0.1)~0)

RRA_result_GSEA_1 <- na.omit(RRA_result_GSEA)

# Heat map generation
pheatmap(RRA_result_GSEA_1, cluster_rows=TRUE , show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col= Name,cellheight = 20, cellwidth= 30)