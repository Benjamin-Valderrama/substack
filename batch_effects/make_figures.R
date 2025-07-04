
# Preparation -------------------------------------------------------------

library(tidyverse)

# metadata together
metadata_together <- read_tsv("batch_effects/inputs/metadata_error_rates_together.tsv")

# ct learn errors together
ct_together <- read_tsv("batch_effects/inputs/ct_error_rates_together.tsv") %>% 
  column_to_rownames("genus")

sort(metadata_together$sample_id) == sort(colnames(ct_together))
metadata_together <- metadata_together[metadata_together$sample_id %in% colnames(ct_together), ]
ct_together <- ct_together[, metadata_together$sample_id]

# metadata separate
metadata_separate <- read_tsv("batch_effects/inputs/metadata_error_rates_separately.tsv")

# ct learn errors separate
ct_separate <- read_tsv("batch_effects/inputs/ct_error_rates_separately.tsv") %>% 
  column_to_rownames("genus")

metadata_separate$sample_id == colnames(ct_separate)
ct_separate <- ct_separate[, metadata_separate$sample_id]


# CLR transformation ------------------------------------------------------
library(vegan)

clr_ct_together <- decostand(x = ct_together, 
                             MARGIN = 2, 
                             method = "clr", 
                             pseudocount = 2/3, )

clr_ct_separate <- decostand(x = ct_separate, 
                             MARGIN = 2, 
                             method = "clr", 
                             pseudocount = 2/3)


# ComBat ------------------------------------------------------------------
library(sva)

mod_together = model.matrix( ~ treatment, data = metadata_together)
ct_together_combat <- ComBat(clr_ct_together, 
                             par.prior = FALSE,
                             batch = metadata_together$batch,
                             mod = mod_together)



mod_separate = model.matrix( ~ treatment, data = metadata_separate)
ct_separate_combat <- ComBat(clr_ct_separate,
                             par.prior = FALSE,
                             batch = metadata_separate$batch,
                             mod = mod_separate)



# visually assess batch effect --------------------------------------------
library(pheatmap)
library(see)

colors <- okabeito_colors(1:4)
names(colors) <- paste0("b", 1:4)
annotation_colors <- list(batch = colors)

together_annotation <- data.frame(batch = metadata_together$batch)
rownames(together_annotation) <- metadata_together$sample_id


plot_together <- pheatmap(mat = clr_ct_together, 
         annotation_col = together_annotation,
         cluster_cols = TRUE, cluster_rows = TRUE, 
         annotation_colors = annotation_colors,
         colorRampPalette(c("dodgerblue3", "#fffff0", "tomato3"))(50), 
         legend = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         #main = "Genera Abundances with Error Rates Learnt Across all Batches"
         )

ggsave(plot = plot_together,
       filename = "batch_effects/outputs/error_rates_together.jpg")


plot_together_combat <- pheatmap(mat = ct_together_combat, 
         annotation_col = together_annotation,
         cluster_cols = TRUE, cluster_rows = TRUE, 
         annotation_colors = annotation_colors,
         colorRampPalette(c("dodgerblue3", "#fffff0", "tomato3"))(50), 
         legend = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         #main = "Genera Abundances with Error Rates Learnt Across all Batches + ComBat"
         )

ggsave(plot = plot_together_combat,
       filename = "batch_effects/outputs/error_rates_together_and_combat.jpg")




separate_annotation <- data.frame(batch = metadata_separate$batch)
rownames(separate_annotation) <- metadata_separate$sample_id


plot_separate <- pheatmap(mat = clr_ct_separate, 
         annotation_col = separate_annotation,
         cluster_cols = TRUE, cluster_rows = TRUE, 
         annotation_colors = annotation_colors,
         colorRampPalette(c("dodgerblue3", "#fffff0", "tomato3"))(50), 
         legend = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         #main = "Genera Abundances with Different Error Rates Learnt for Each Batch"
         )

ggsave(plot = plot_separate,
       filename = "batch_effects/outputs/error_rates_separate.jpg")



plot_separate_combat <- pheatmap(mat = ct_separate_combat, 
         annotation_col = separate_annotation,
         cluster_cols = TRUE, cluster_rows = TRUE, 
         annotation_colors = annotation_colors,
         colorRampPalette(c("dodgerblue3", "#fffff0", "tomato3"))(50), 
         legend = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         #main = "Genera Abundances with Different Error Rates Learnt for Each Batch + ComBat"
         )

ggsave(plot = plot_separate_combat,
       filename = "batch_effects/outputs/error_rates_separate_and_combat.jpg")