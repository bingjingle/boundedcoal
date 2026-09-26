# Invoke with Rscript; archived input paths may be overridden with an environment variable.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_file <- gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)
script_dir <- dirname(normalizePath(script_file))
data_dir <- Sys.getenv("BOUNDEDCOAL_CELL_LINEAGE_DATA_DIR", script_dir)
output_dir <- Sys.getenv("BOUNDEDCOAL_CELL_LINEAGE_OUTPUT_DIR", file.path(script_dir, "../../outputs/cell_lineage"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load R packages
library(ape)
library(phangorn)
library(tidyverse)
library(ggdendro)
library(dendextend)

options(digits = 22)
# =====================================================================
# Generating the wider cell-by-59EditSites file from cell-by-5EditSites
# =====================================================================
edit_table_by_5 = read.csv(file.path(data_dir, "Supplementary_File_2_DataTableMOI19.csv"), stringsAsFactors = F, header = T, na.strings=c("","NA"))

### Filtering TargetBC set based on the number of retained cell (13 TargetBC, 3257 cells)
TargetBC_freq <- data.frame(table(edit_table_by_5$TargetBC)) %>% arrange(desc(Freq))
retained_cells <- rep(0,17)
for (ii in 1:17){
  current_TargetBC_set <- as.character(TargetBC_freq$Var1)[1:ii]
  Cell_filter_by_TargetBC <- select(edit_table_by_5, c('Cell','TargetBC')) %>%
    filter(TargetBC %in% current_TargetBC_set) %>%
    mutate(count = 1) %>%
    pivot_wider(names_from = TargetBC, values_from = 'count', values_fill = 0)
  retained_cells[ii] <- sum(rowSums(Cell_filter_by_TargetBC[,2:(ii+1)]) == ii)
}
retained_cells <- data.frame(retained_cells)
retained_cells$numTargetBC <- 1:17
retained_cells$retained_cells <- as.numeric(as.character(retained_cells$retained_cells))
targets_to_use <- as.character(TargetBC_freq$Var1[1:13])


### Generating the cell_list for 3,257 cells
Cell_filter_13set <- select(edit_table_by_5, c('Cell','TargetBC')) %>%
  filter(TargetBC %in% targets_to_use) %>%
  mutate(count = 1) %>%
  pivot_wider(id_cols = 'Cell',names_from = TargetBC, values_from = 'count', values_fill = 0)
Cell_filter_13set$sums <- rowSums(Cell_filter_13set[,2:14])
Cell_filter_13set <- filter(Cell_filter_13set, sums == 13)
cell_list_3257 <- Cell_filter_13set$Cell
edit_table_3257 <- filter(edit_table_by_5, TargetBC %in% targets_to_use & Cell %in% cell_list_3257)

# Marking unedited sites with 'None', and non-existing sites (contracted to 4xTAPE or 2xTAPE) with NA
edit_table_3257[is.na(edit_table_3257)] <- 'None'
edit_table_3257[edit_table_3257$TargetBC == 'TGGACGAC',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TTTCGTGA',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TGGTTTTG',7] <- NA
edit_table_3257[edit_table_3257$TargetBC == 'TTCACGTA',5:7] <- NA


# edit_cell_table_65 = Ordered cell-by-65EditSites table, including non-existing sites before contraction
edit_cell_table_65 <- select(edit_table_3257, -nUMI) %>%
  pivot_longer(cols = c('Site1','Site2','Site3','Site4','Site5'), names_to = 'Sites', values_to ='Insert') %>%
  pivot_wider(id_cols = Cell, names_from = c(TargetBC,Sites), names_sep = ".", values_from = Insert)
edit_cell_table_65 <- arrange(edit_cell_table_65,Cell) %>%
  select(order(colnames(edit_cell_table_65)))


sub_edit65 <- as.matrix(select(edit_cell_table_65,-Cell))
sub_edit65[is.na(sub_edit65)] <- 'None'
rownames(sub_edit65) <- edit_cell_table_65$Cell

sub_edit59 <- as.matrix(select(edit_cell_table_65,-c('Cell','TGGACGAC.Site5','TTTCGTGA.Site5','TGGTTTTG.Site5',
                                                     'TTCACGTA.Site3','TTCACGTA.Site4','TTCACGTA.Site5')))
rownames(sub_edit59) <- edit_cell_table_65$Cell
cell_list <- edit_cell_table_65$Cell


# =====================================================================
# Generating the phylogenetic tree based on edits
# =====================================================================

# shared_edit_matrix = Counting all shared edits per cell-pair, consistent with the sequential editing on DNA Tape



# Function for calculating shared_edit_matrix
fun_shared_edit_matrix <- function(x) {
  #if (1 == 1){
  sub_edit65 <- x
  sub_edit65[sub_edit65 == 'None'] <- 1:filter(as.data.frame(table(sub_edit65)), sub_edit65 == 'None')$Freq
  ncell <- nrow(sub_edit65)
  cell_list <- sort(rownames(sub_edit65))
  shared_edit_matrix <- matrix(0, ncell,ncell)
  colnames(shared_edit_matrix) <- cell_list
  rownames(shared_edit_matrix) <- cell_list
  for (ii in 1:(ncell)){
    cell1 <- cell_list[ii]
    for (jj in (ii):ncell){
      cell2 <- cell_list[jj]
      for (kk in seq(0,(dim(sub_edit65)[2]-5),5)){
        if (sub_edit65[cell1,(kk+1)] == sub_edit65[cell2,(kk+1)]){
          shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
          if (sub_edit65[cell1,(kk+2)] == sub_edit65[cell2,(kk+2)]){
            shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
            if (sub_edit65[cell1,(kk+3)] == sub_edit65[cell2,(kk+3)]){
              shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
              if (sub_edit65[cell1,(kk+4)] == sub_edit65[cell2,(kk+4)]){
                shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                if (sub_edit65[cell1,(kk+5)] == sub_edit65[cell2,(kk+5)]){
                  shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                }
              }
            }
          }
        }
      }
      shared_edit_matrix[jj,ii] <- shared_edit_matrix[ii,jj]
    }
  }
  return(shared_edit_matrix)
}

#shared_edit_matrix <- fun_shared_edit_matrix(sub_edit65) #15-30 min on computer; once done, saved and loaded for the future use
shared_edit_matrix <- read.csv(file.path(data_dir, "shared_edit_matrix_3257.csv"), stringsAsFactors = F, header = T)
shared_edit_matrix <- as.matrix(shared_edit_matrix)
diag(shared_edit_matrix) <- 59

distance_matrix <- 59 - shared_edit_matrix # Phylogenetic distance caludated as (# of possible sites - # of shared sites)
distance_matrix <- as.matrix(distance_matrix)
tree <- as.phylo(hclust(as.dist(distance_matrix), "average")) # tree built using UPGMA
tree
tree_height_calc <- function(tree) {
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]

  }
  return(sum_path)

}
tree_height <- tree_height_calc(tree)
tree$edge.length <- tree$edge.length * (25.0/tree_height)
tree$edge.length <- tree$edge.length - 0.001
bt <- sort(branching.times(tree))

min_gap <- 1e-4
bt_adj <- bt

for (i in 2:length(bt_adj)) {
  if (bt_adj[i] - bt_adj[i - 1] < min_gap) {
    bt_adj[i] <- bt_adj[i - 1] + min_gap
  }
}

save(bt_adj, file=file.path(output_dir, "UGPMA3257new.rda"))











########################################################################################
tree <- read.tree(file = file.path(data_dir, "UPGMAtree_100.txt"))
bt <- sort(branching.times(tree))
min_gap <- 1e-4
bt_adj <- bt
for (i in 2:length(bt_adj)) {
  if (bt_adj[i] - bt_adj[i - 1] < min_gap) {
    bt_adj[i] <- bt_adj[i - 1] + min_gap
  }
}
save(bt_adj, file=file.path(output_dir, "UGPMA100new.rda"))
