library(xtable)

tt <- names(tt) %>%
    lapply(\(x){
        tt[[x]] %>%
            dplyr::filter(chr %in% primary_chrs) %>%
            left_join(as_tibble(gene_dar[[x]])[,c("gene_id", "dar")])
    }) %>%
    set_names(names(tt))

latex_table <- function(tt){
    tt %>%
        dplyr::filter(DE) %>%
        mutate(
            PValue = formatC(PValue, digits = 2, format = "e"),
            FDR = formatC(FDR, digits = 2, format = "e")
        ) %>%
        dplyr::select(
            Gene = gene_name, ID = gene_id, Chromosome = chr, logFC, p = PValue, FDR, DAR = dar
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

latex_table_post <- function(tt){
    tt %>%
        dplyr::filter(darDE) %>%
        mutate(
            PValue = formatC(darP, digits = 2, format = "e"),
            FDR = formatC(darFDR, digits = 2, format = "e")
        ) %>%
        dplyr::select(
            Gene = gene_name, ID = gene_id, Chromosome = chr, logFC, p = PValue, FDR, DAR = dar
        ) %>%
        xtable() %>%
        print(include.rownames = FALSE)
}

latex_table_post(tt_dar$pk_tu)
latex_table_post(tt_dar$pk_het)
# latex_table_post(tt$tu_het)
latex_table_post(tt_dar$eofad_wt)
latex_table_post(tt_dar$fai_wt)
latex_table_post(tt_dar$eofad_fai)
