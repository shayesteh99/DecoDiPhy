library(jsonlite)
library(ggtree)
library(stringr)
library(dplyr)
library(tidyverse)
library(treeio)
library(viridisLite)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
mode = args[2]
output_file <- args[3]

# input_file <- "../Example/output/output.jplace"
# mode = "weighted"
# output_file <- "../Example/output/DecoDiPhy_results.pdf"

jp <- fromJSON(input_file, simplifyDataFrame = FALSE)

#combine reads
edge_weights <- list()

for (p in jp$placements) {
  total = nrow(p$p)
  for (i in seq_len(total)) {
    edge <- p$p[i, 1]   
    if (mode == "unweighted") {
      weight <- 1/total
    } else {
      weight <- p$p[i, 2]
    }
    
    edge_weights[[as.character(edge)]] <-
      (edge_weights[[as.character(edge)]] %||% 0) + weight
  }
}

edge_df <- data.frame(
  id = as.integer(names(edge_weights)),
  abundance = unlist(edge_weights)
)

read_count = sum(edge_df$abundance)
edge_df$abundance <- edge_df$abundance/read_count
edge_df$id <- edge_df$id

#read tree
nodes <- str_match_all(jp$tree, "([A-Za-z0-9_]+):[0-9\\.Ee-]+\\{([0-9]+)\\}")[[1]]
nodes_df <- data.frame(
  label = nodes[,2],
  id = as.integer(nodes[,3])
)

tree_phy <- read.tree(text = jp$tree)
p <- ggtree(tree_phy)
tree_df <- p$data

tree_df <- merge(tree_df, nodes_df)
edge_df <- merge(edge_df, tree_df, by = "id")

ggtree(tree_phy) %<+% edge_df +
  geom_tree(aes(color = abundance, linewidth=is.na(abundance))) +
  scale_linewidth_manual(values = c(1.5, 0.5), guide="none")+
  scale_color_gradientn(
    colours = viridis(256, option = "C"),
    trans = "log10",
    na.value = "black",
    limits = c(1e-9, 3e-1),
    name = "abundance"
  ) +
  geom_text2(aes(label = label,subset = startsWith(label, "G")), 
             color="black", nudge_x = 0.05, size = 2)+
  theme_tree()

ggsave(output_file, width= 4.5, height = 11)
