## Download the seurat obj from 10.5281/zenodo.15747708

library(Seurat)
library(ggplot2)
library(ggforce) 

# List of genes of interest for plotting
genes_of_interest <- "Cd3e" # You can add other genes, e.g. c("Cd3e", "Cd4", "Foxp3")
violin_plots <- list()

# Define the y-axis limits based on the global maximum value across the genes of interest
y_limits <- c(
  0,
  max(
    sapply(
      genes_of_interest,
      function(gene) max(FetchData(seurat.integrated, vars = gene), na.rm = TRUE)
    )
  )
)

plots_dir <- file.path(path, "Plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

for (gene in genes_of_interest) {
  gene_data <- FetchData(seurat.integrated, vars = c("condition", gene))
  gene_data$cluster <- Idents(seurat.integrated)
  
  # Generate the violin plot
  violin_plot <- ggplot(gene_data, aes(x = condition, y = !!sym(gene), fill = condition)) +
    geom_sina(aes(color = condition), size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3, color = "black", alpha = 0.6) +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_fill_manual(values = c("notch3" = "red", "wt" = "black")) +
    scale_color_manual(values = c("notch3" = "red", "wt" = "black")) +
    coord_cartesian(ylim = y_limits) +
    theme_classic(base_family = "Arial") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 18, face =  "bold"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      strip.background = element_rect(color = "black", fill = "grey90"),
      strip.text = element_text(size = 20, face = "bold")
    ) +
    labs(title = gene, x = "Condition", y = "Expression Level") +
    scale_x_discrete(limits = c("wt", "notch3"))
  
  violin_plots[[gene]] <- violin_plot
  
  # Salva il plot
  ggsave(filename = file.path(plots_dir, paste0("violin_", gene, ".pdf")),
         plot = violin_plot, width = 10, height = 8)
}
