library(patchwork)
library(ggside)

theme_update(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    title = element_markdown()
)

pal_paired <- brewer.pal(8, "Paired")

seq_lengths <- seqlengths(genes) %>%
    .[primary_chrs]

tt <- readRDS(here("Rdata/tt.Rds"))
gene_dar <- readRDS(here("Rdata/gene_dar_1e6.Rds"))

psen1_midpoint <- genes %>%
    as_tibble() %>%
    dplyr::filter(gene_name == "psen1") %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["17"]) %>%
    pull(midpoint)

## D/D vs A/B

de_midpoint <- tt$pk_tu %>%
    dplyr::filter(DE, chr == 14) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["14"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$pk_tu) %>%
    left_join(as_tibble(tt$pk_tu)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 14) %>%
    pull(dar)
dar_1e6$pk_tu %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_xsidevline(xintercept = de_midpoint, linetype = "solid", colour = "red") +
    geom_ysidehline(yintercept = de_dar, linetype = "solid", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *PK<sup>D</sup>/PK<sup>D</sup>* vs. *TU<sup>A</sup>/TU<sup>B</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        # ggside.panel.scale.x = 0.025,
        # ggside.panel.scale.y = 0.01875
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr14_pktu.png",
    width = 8, height = 4
)

## D/D vs B/D & A/B vs B /D

de_midpoint <- tt$pk_het %>%
    dplyr::filter(DE, chr == 14) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["14"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$pk_het) %>%
    left_join(as_tibble(tt$pk_het)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 14) %>%
    pull(dar)
a <- dar_1e6$pk_het %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_xsidevline(xintercept = de_midpoint, linetype = "solid", colour = "red") +
    geom_ysidehline(yintercept = de_dar, linetype = "solid", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *PK<sup>D</sup>/PK<sup>D</sup>* vs. *TU<sup>B</sup>/PK<sup>D</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "right",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
de_midpoint <- tt$tu_het %>%
    dplyr::filter(DE, chr == 14) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["14"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$tu_het) %>%
    left_join(as_tibble(tt$tu_het)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 14) %>%
    pull(dar)
b <- dar_1e6$tu_het %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_xsidevline(xintercept = 0.5, linetype = "blank", colour = "red") +
    geom_ysidehline(yintercept = 0.5, linetype = "blank", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *TU<sup>A</sup>/TU<sup>B</sup>* vs. *TU<sup>B</sup>/PK<sup>D</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "right",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
a + b +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
        legend.position = "bottom",
        plot.margin = unit(c(.2, .5, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    )
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr14_hets.png",
    width = 8, height = 8
)

## T4/+ vs +/+ & W2/+ vs +/+

de_midpoint <- tt$eofad_wt %>%
    dplyr::filter(DE, chr == 17) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["17"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$eofad_wt) %>%
    left_join(as_tibble(tt$eofad_fai)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 17) %>%
    pull(dar)
c <- dar_1e6$eofad_wt %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_vline(xintercept = psen1_midpoint, linetype = "dashed") +
    geom_ysidehline(yintercept = de_dar, linetype = "solid", colour = "red") +
    geom_text(x = psen1_midpoint - .05, y = 0.95, label = "psen1", fontface = "italic") +
    geom_xsidevline(xintercept = de_midpoint, linetype = "solid", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>T428del</sup>/PK<sup>+</sup>* vs. *PK<sup>+</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "right",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
de_midpoint <- tt$fai_wt %>%
    dplyr::filter(DE, chr == 17) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["17"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$fai_wt) %>%
    left_join(as_tibble(tt$fai_wt)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 17) %>%
    pull(dar)
d <- dar_1e6$fai_wt %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_vline(xintercept = psen1_midpoint, linetype = "dashed") +
    geom_text(x = psen1_midpoint - .05, y = 0.95, label = "psen1", fontface = "italic") +
    geom_xsidevline(xintercept = de_midpoint, linetype = "solid", colour = "red") +
    geom_ysidehline(yintercept = de_dar, linetype = "solid", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>W233fs</sup>/PK<sup>+</sup>* vs. *PK<sup>+</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "right",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
c + d +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
        legend.position = "bottom",
        plot.margin = unit(c(.2, .5, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    )
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr17_eofadwt_faiwt.png",
    width = 8, height = 8
)

## T4/PK+ vs W2/PK+

de_midpoint <- tt$eofad_fai %>%
    dplyr::filter(DE, chr == 17) %>%
    mutate(midpoint = ((start + end) / 2) / seq_lengths["17"]) %>%
    pull(midpoint)
de_dar <- as_tibble(gene_dar$eofad_fai) %>%
    left_join(as_tibble(tt$eofad_fai)[,c("gene_id", "chr", "DE")]) %>%
    dplyr::filter(DE, chr == 17) %>%
    pull(dar)
dar_1e6$eofad_fai %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    geom_vline(xintercept = psen1_midpoint, linetype = "dotted") +
    geom_text(x = psen1_midpoint - .05, y = 0.95, label = "psen1", fontface = "italic") +
    geom_xsidevline(xintercept = de_midpoint, linetype = "solid", colour = "red") +
    geom_ysidehline(yintercept = de_dar, linetype = "solid", colour = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = pal_paired[6], "ln_true" = pal_paired[6],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>T428del</sup>/PK<sup>+</sup>* vs. *TU<sup>W233fs</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr17_eofadfai.png",
    width = 8, height = 4
)




####
## SECONDARY LOCUS
####

dar_1e6$fai_wt %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>W233fs</sup>/PK<sup>+</sup>* vs. *PK<sup>+</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr14_faiwt.png",
    width = 8, height = 4
)


e <- dar_1e6$eofad_wt %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>T428del</sup>/PK<sup>+</sup>* vs. *PK<sup>+</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
f <- dar_1e6$eofad_fai %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "14", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "14", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "*psen1* *TU<sup>T428del</sup>/PK<sup>+</sup>* vs. *TU<sup>W233fs</sup>/PK<sup>+</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 14"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
e + f +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
        legend.position = "bottom",
        plot.margin = unit(c(.2, .5, .2, .2), "cm")
    )
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr14_eofadwt_eofadfai.png",
    width = 8, height = 8
)


g <- dar_1e6$pk_tu %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *PK<sup>D</sup>/PK<sup>D</sup>* vs. *TU<sup>A</sup>/TU<sup>B</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
h <- dar_1e6$tu_het %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *TU<sup>A</sup>/TU<sup>B</sup>* vs. *TU<sup>A</sup>/PK<sup>D</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
i <- dar_1e6$pk_het %>%
    as_tibble() %>%
    split(.$seqnames) %>%
    lapply(\(x){
        chr <- unique(x$seqnames)
        chr_length <- seq_lengths[as.character(chr)]
        x %>%
            mutate(rel_position = start /chr_length)
    }) %>%
    bind_rows() %>%
    mutate(
        point_group = ifelse(seqnames == "17", "pt_true", "pt_false"),
        point_group = fct_relevel(as.character(point_group), "pt_false"),
        line_group = ifelse(seqnames == "17", "ln_true", "ln_false"),
        line_group = fct_relevel(as.character(line_group), "ln_false"),
    ) %>%
    dplyr::arrange(point_group) %>%
    ggplot(aes(rel_position, dar_region)) +
    geom_point(aes(colour = point_group), size = 0.5, show.legend = FALSE) +
    geom_line(aes(colour = point_group), data = . %>% dplyr::filter(point_group == "pt_true"), show.legend = TRUE) +
    # geom_smooth(aes(colour = line_group), se = FALSE) +
    geom_smooth(aes(colour = line_group), data = . %>% dplyr::filter(line_group == "ln_false"), se = FALSE, linewidth = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.51)) +
    scale_colour_manual(
        values = c(
            "pt_true" = "grey40", "ln_true" = pal_paired[8],
            "pt_false" = pal_paired[1], "ln_false" = pal_paired[2]
        ),
        labels = c("pt_true" = "True", "ln_false" = "False"),
        breaks = c("pt_true", "ln_false")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        title = "chr14 *PK<sup>D</sup>/PK<sup>D</sup>* vs. *TU<sup>A</sup>/PK<sup>D</sup>*",
        x = "Relative Chromosomal Position",
        y = "DAR (region)",
        colour = "Chromosome 17"
    ) +
    theme_pubclean() +
    theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(.2, 1, .2, .2), "cm"),
        ggside.panel.scale.x = 0.075,
        ggside.panel.scale.y = 0.03
    ) +
    ggside(x.pos = "bottom", y.pos = "left")
g + h + i +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
        legend.position = "bottom",
        plot.margin = unit(c(.2, .5, .2, .2), "cm")
    )
ggsave(
    "~/phd/publications/pktu_manuscript/fig/dar_chr17_pktu_tuhet_pkhet.png",
    width = 8, height = 12
)

