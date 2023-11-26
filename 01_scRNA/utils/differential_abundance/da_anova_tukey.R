# Test difference in proportion of cell types by one-way ANOVA with Tukey’s multiple comparison

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(ggpubr)
})

# set paths
path_main <- '/Users/ieo4874/Desktop/working/SCMseq/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/all_samples/06_lineage_analysis/differential_abundance/')

# settings
load(paste0(path_data, "settings.RData"))

# define function
run_anova_tukey <- function(data, seurat) {
  tab_prop <- left_join(
    data, 
    distinct(seurat@meta.data[,c("sample", "cohort")]), 
    by = "sample"
  ) %>%
    relocate("cohort", .after = "sample") %>%
    pivot_longer(cols = 4:ncol(.), names_to = "lineage", values_to = "count") %>%
    mutate(proportion = 100*(count / total_cell_count),
           lineage = str_replace_all(lineage, "_", " "))
  
  for (l in unique(tab_prop$lineage)) {
    x <- tab_prop %>%
      filter(lineage == l) 
    
    # test equality of variances
    res.aov <- aov(proportion ~ cohort, data = x)
    levene_test <- car::leveneTest(proportion ~ cohort, data = x)
    
    # homogeneity of variances in different groups can be assumed if p > 0.05
    if (levene_test[[3]][1] > 0.05) {
      
      # run Tukey
      p = ggplot(x, aes(cohort, proportion, fill = cohort)) +
        geom_jitter(size = 2.5, width = 0.3, height = 0,
                    shape = 21, stroke = 0.2,
                    colour = "black") +
        geom_point(data = x %>% 
                     group_by(cohort) %>% 
                     summarise(mean = mean(proportion)), 
                   mapping = aes(x = cohort, y = mean), 
                   size = 10, color = 'black', shape = '_') +
        theme_bw() +
        scale_y_continuous(expand = c(0.1,-0.1)) +
        ylab("% cells") +
        ggtitle(l) +
        scale_fill_manual(values = colors$cohort) +
        stat_pwc(method = "tukey_hsd",
                 bracket.shorten = 0.3,
                 tip.length = 0,
                 vjust = 0,
                 label.size = 2.8,
                 label = "p.adj.signif",
                 symnum.args = symnum.args,
                 hide.ns = T
        ) +
        theme(panel.grid = element_blank(),
              aspect.ratio = 1,
              legend.position = "none", 
              axis.title.x = element_blank(),
              axis.text = element_text(color = "black"),
              axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
              axis.title.y = element_text(size = 10),
              plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2, size = 12),
              plot.margin = unit(c(1,1,0,0), "cm")
        )
      ggsave(paste0(path_results, "plots/", l, "_anova_tukey.png"), p, width = 4, height = 4)
      
    } else {
      
      message("Anova not applicable for lineage: ", l)
      next 
      
    }
  }
}
