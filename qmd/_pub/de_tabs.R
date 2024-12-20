library(xtable)

latex_table <- function(tt){
    tt %>%
        dplyr::filter(DE) %>%
        mutate(
            PValue = formatC(PValue, digits = 2, format = "e"),
            FDR = formatC(FDR, digits = 2, format = "e")
        ) %>%
        dplyr::select(
            Gene = gene_name, ID = gene_id, Chromosome = chr, logFC, p = PValue, FDR
        ) %>%
        xtable() %>%
        print(include.rownames = FALSE)
}

latex_table(tt$pk_tu)
latex_table(tt$pk_het)
# latex_table(tt$tu_het)
latex_table(tt$eofad_wt)
latex_table(tt$fai_wt)
latex_table(tt$eofad_fai)
