#FGSEA Analysis
library('tidyverse')
library( "biomaRt" )
library('fgsea')
library('qusage')
Pathways <- gmtPathways('/Users/ikjotsidhu/Downloads/c2.cp.v6.2.symbols.gmt')

# Getting gene symbols
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )


# Gene symbols for HD_DE genes
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = HD_DE$X1,
                  mart = ensembl )

idx <- match( HD_DE$X1, genemap$ensembl_gene_id )

HD_DE$entrez <- genemap$entrezgene[ idx ]
HD_DE$SYMBOL <- genemap$hgnc_symbol[ idx ]

# Gene symbols for PD_DE genes
genemap1 <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = PD_DE$X1,
                  mart = ensembl )

idx <- match( PD_DE$X1, genemap1$ensembl_gene_id )

PD_DE$entrez <- genemap1$entrezgene[ idx ]
PD_DE$SYMBOL <- genemap1$hgnc_symbol[ idx ]

# Preparing HD data for fgsea

signif_res_ranks_HD <- HD_DE %>% 
  select(symbol, cnts.beta) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(cnts.beta=mean(cnts.beta))

signif_res_ranks_HD <- signif_res_ranks_HD[!(is.na(signif_res_ranks_HD$symbol) | signif_res_ranks_HD$symbol==""), ]

ranks_HD <- deframe(signif_res_ranks_HD)

# Preparing PD data for fgsea

signif_res_ranks_PD <- PD_DE %>% 
  select(symbol, cnts.beta) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(cnts.beta=mean(cnts.beta))

signif_res_ranks_PD <- signif_res_ranks_PD[!(is.na(signif_res_ranks_PD$symbol) | signif_res_ranks_PD$symbol==""), ]

ranks_PD <- deframe(signif_res_ranks_PD)

# Preparing CTES vs RHIN data for fgsea
signif_res_ranks_CTESvRHIN <- CTES_RHIN_DE %>% 
  select(gene_name,Status__CTES__log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(Status__CTES__log2FoldChange=mean(Status__CTES__log2FoldChange))

#signif_res_ranks_PD <- signif_res_ranks_PD[!(is.na(signif_res_ranks_PD$symbol) | signif_res_ranks_PD$symbol==""), ]

ranks_CTESvRHIN <- deframe(signif_res_ranks_CTESvRHIN)

# Preparing CTEM vs RHIN data for fgsea
signif_res_ranks_CTEMvRHIN <- CTEM_RHIN_DE %>% 
  select(gene_name,Status__CTEM__log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(Status__CTEM__log2FoldChange=mean(Status__CTEM__log2FoldChange))

#signif_res_ranks_PD <- signif_res_ranks_PD[!(is.na(signif_res_ranks_PD$symbol) | signif_res_ranks_PD$symbol==""), ]

ranks_CTEMvRHIN <- deframe(signif_res_ranks_CTEMvRHIN)



#Run fgsea for HD
fgseaRes_HD <- fgsea(pathways=Pathways, stats=ranks_HD, nperm=1000)
fgseaResTidy_HD <- fgseaRes_HD %>%
  as_tibble() %>%
  arrange(pval)

#Run fgsea for PD
fgseaRes_PD <- fgsea(pathways=Pathways, stats=ranks_PD, nperm=1000)
fgseaResTidy_PD <- fgseaRes_PD %>%
  as_tibble() %>%
  arrange(pval)

#Run fgsea for CTESvsRHIN 
fgseaRes_CTESvsRHIN <- fgsea(pathways=Pathways, stats=ranks_CTESvRHIN, nperm=1000)
fgseaResTidy_CTESvsRHIN <- fgseaRes_CTESvsRHIN %>%
  as_tibble() %>%
  arrange(desc(NES))

