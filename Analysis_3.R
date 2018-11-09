# Import DE results from HD,PD and CTES vs RHIN

#Sort CTESvsRHIN and CTEMvsRHIN statistics by Status__CTES__pvalue/Status__CTEM__pvalue and add rank columns

CTES_RHIN_DE <- read_csv("CTES_RHIN_all_filthi_de_annot.csv") %>% filter(!is.na(gene_name))
CTES_RHIN_DE <-  arrange(Status__CTES__pvalue) %>% mutate(CTE_RHIN_Rank=row_number())

CTEM_RHIN_DE <- read_csv("CTEM_RHIN_all_filthi_de_annot.csv") %>% filter(!is.na(gene_name))
CTEM_RHIN_DE <-  arrange(Status__CTEM__pvalue) %>% mutate(CTE_RHIN_Rank=row_number())

# Sort the HD and PD DE statistics by cnts.p and add rank columns
HD_DE <- read_csv("HD_DE.csv") %>% filter(!is.na(symbol)) %>% arrange(cnts.p) %>% mutate(HD_rank=row_number())
PD_DE <- read_csv("PD_DE.csv") %>% filter(!is.na(symbol)) %>% arrange(cnts.p) %>% mutate(PD_rank=row_number())

# Perform RRA

# Get the sorted genes
HD_DE_genes_sorted <- (HD_DE$symbol)
PD_DE_genes_sorted <- (PD_DE$symbol)
CTEM_DE_genes_sorted <- (CTEM_RHIN_DE$gene_name)
CTES_DE_genes_sorted <- (CTES_RHIN_DE$gene_name)

#Pass the sorted genes to a single list for RRA

DE_sorted_list <- list(HD_DE_genes_sorted,PD_DE_genes_sorted,CTEM_DE_genes_sorted,CTES_DE_genes_sorted)

ranked_DE_genes <- aggregateRanks(DE_sorted_list)

ranked_DE_genes <- remove_rownames(ranked_DE_genes)

RRA_result_DE_genes<- ranked_DE_genes %>% arrange(Score)


# ADDING COLUMNS
#1. the rank
#2. cnts.beta (HD,PD) or Status__CTES__log2FoldChange/Status__CTEM__log2FoldChange (CTE)
#3. cnts.padj (HD,PD) or Status__CTES__padj/Status__CTEM__padj (CTE)
#of each gene from each of HD, PD, and CTE to the dataframe

# Matching cnts.beta from HD
RRA_result_DE_genes$cnts.beta_HD <- HD_DE$cnts.beta[match( RRA_result_DE_genes$Name,HD_DE$symbol)]

# Matching cnts.beta from PD
RRA_result_DE_genes$cnts.beta_PD <- PD_DE$cnts.beta[match( RRA_result_DE_genes$Name,PD_DE$symbol)]

# Matching Status__CTES__log2FoldChange from CTES
RRA_result_DE_genes$Status__CTES__log2FoldChange <- CTES_RHIN_DE$Status__CTES__log2FoldChange[match(RRA_result_DE_genes$Name, CTES_RHIN_DE$gene_name)]

# Matching Status__CTEM__log2FoldChange from CTEM
RRA_result_DE_genes$Status__CTEM__log2FoldChange <- CTEM_RHIN_DE$Status__CTEM__log2FoldChange[match(RRA_result_DE_genes$Name, CTEM_RHIN_DE$gene_name)]

# Matching padj from HD
RRA_result_DE_genes$padj_HD <- HD_DE$cnts.padj[match( RRA_result_DE_genes$Name,HD_DE$symbol)]

# Matching padj from PD
RRA_result_DE_genes$padj_PD <- PD_DE$cnts.padj[match( RRA_result_DE_genes$Name,PD_DE$symbol)]

# Matching Status__CTES__padj from CTES
RRA_result_DE_genes$Status__CTES__padj <- CTES_RHIN_DE$Status__CTES__padj[match(RRA_result_DE_genes$Name, CTES_RHIN_DE$gene_name)]

# Matching Status__CTES__padj from CTES
RRA_result_DE_genes$Status__CTEM__padj <- CTEM_RHIN_DE$Status__CTEM__padj[match(RRA_result_DE_genes$Name, CTEM_RHIN_DE$gene_name)]

# Matching rank from HD
RRA_result_DE_genes$HD_rank <- HD_DE$HD_rank[match( RRA_result_DE_genes$Name,HD_DE$symbol)]

# Matching padj from PD
RRA_result_DE_genes$PD_rank <- PD_DE$PD_rank[match( RRA_result_DE_genes$Name,PD_DE$symbol)]

# Matching Status__CTES__padj from CTES
RRA_result_DE_genes$CTES_rank <- CTES_RHIN_DE$CTE_RHIN_Rank[match(RRA_result_DE_genes$Name, CTES_RHIN_DE$gene_name)]

# Matching Status__CTES__padj from CTES
RRA_result_DE_genes$CTEM_rank <- CTEM_RHIN_DE$CTE_RHIN_Rank[match(RRA_result_DE_genes$Name, CTEM_RHIN_DE$gene_name)]

RRA_result_DE_genes <- na.omit(RRA_result_DE_genes)
write_csv(RRA_result_DE_genes,"RRA_result_DE.csv")