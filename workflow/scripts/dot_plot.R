# Check if running as Snakemake script or command line
if (exists("snakemake")) {
  # Running as Snakemake script
  input_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]
} else {
  # Running from command line
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript dot_plot.R <input_file> <output_file>")
  }
  input_file <- args[1]
  output_file <- args[2]
}

library(ggplot2)
library(ggpubr)

# Read the data
df <- read.table(input_file, sep = '\t', header = TRUE)

# Create dot plot for contigs
dp <- ggplot(df, aes(x=who, y=gc_content, fill=who)) +
    geom_dotplot(binaxis='y', stackdir='up',
        stackratio=0.025, dotsize=0.5) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "bottom") +
    labs(title = "GC Content Distribution by Organism Type",
        x = "Organism Type", 
        y = "GC Content (%)",
        fill = "Type")

# Create histogram for overall GC distribution
h <- ggplot(df, aes(x=gc_content)) +
    geom_histogram(binwidth=2, alpha=0.7, color="black", 
        fill="lightblue") +
    geom_density(alpha=0.3, fill="#FF6666") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(title="Overall GC Content Distribution",
        x = "GC Content (%)", 
        y = "Frequency")

# Combine plots
combined_plot <- ggarrange(h, dp, ncol = 1, nrow = 2)

# Save the plot
ggsave(output_file, plot = combined_plot, width = 10, height = 8, dpi = 300)
