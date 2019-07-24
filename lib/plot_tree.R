#' plot_trehalose_tree
#' Plot the trehalose project tree with annotations for ribotypes, trehalose
#'   utilization variants, and severe infection outcome. Saves plot as a PDF
#'   in the ../figures/ dir. 
#' @param tree_path Character. Path to location of tree file. 
#' @param metadata_path Character. Path to location of metadata file. 
#' @param show_labels Logical. Whether or not to plot tree with isolate names.
#'
#' @noRd
plot_trehalose_tree <- function(tree_path, metadata_path, show_labels = TRUE){
  tree <- read.tree(tree_path)
  metadata <- read_tsv(metadata_path, 
                       col_names = TRUE)
  # GROUPS
  ribotype_014 <- 
    metadata %>% filter(Ribotype == "014-020") %>% select(ID) %>% unlist
  ribotype_015 <- 
    metadata %>% filter(Ribotype == "015") %>% select(ID) %>% unlist
  ribotype_017 <-
    metadata %>% filter(Ribotype == "017") %>% select(ID) %>% unlist
  ribotype_027 <-
    metadata %>% filter(Ribotype == "027") %>% select(ID) %>% unlist
  ribotype_053 <- 
    metadata %>% filter(Ribotype == "053-163") %>% select(ID) %>% unlist
  ribotype_078 <- 
    metadata %>% filter(Ribotype == "078-126") %>% select(ID) %>% unlist
  ribotype_other <- 
    metadata %>% filter(Ribotype == "other") %>% select(ID) %>% unlist
  cases <- metadata %>% filter(Severe_Outcome == 1) %>% select("ID") %>% unlist
  four_gene_insertion <- 
    metadata %>% filter(four_gene_insertion == 1) %>% select("ID") %>% unlist
  L172I <- metadata %>% filter(Leu172Ile == 1) %>% select("ID") %>% unlist
  C171S <- metadata %>% filter(Cys171Ser == 1) %>% select("ID") %>% unlist
  
  # LABEL OFFSET
  ribotype_offset <- 0.001
  cases_offset <- .002
  L172I_offset <- .003
  C171S_offset <- .004
  four_offset <- .004
  
  if (show_labels) {
    base_offset <- 0.001
    ribotype_offset <- base_offset + .003
    cases_offset <- base_offset + .004
    L172I_offset <- base_offset + .005
    C171S_offset <- base_offset + .006
    four_offset <-  base_offset + .006
  }
  
  # SHAPES
  circle <- 19
  triangle <- 17
  empty_circle <- 1
  
  # SIZES 
  dot_size <- .55
  
  # COLORS 
  cases_color <- "red"
  L172I_color <- "deepskyblue"
  C171S_color <- "coral"
  four_color <- "black"
  ribotype_027_color <- "blue"
  ribotype_other_color <- "grey"
  ribotype_014_color <- "deeppink"
  ribotype_015_color <- "orange"
  ribotype_017_color <- "darkorchid"
  ribotype_053_color <- "aquamarine"
  ribotype_078_color <- "darkgreen"
  
  # LEGEND
  legend_colors <- c(L172I_color,
                     C171S_color, 
                     four_color, 
                     cases_color, 
                     ribotype_014_color, 
                     ribotype_015_color, 
                     ribotype_017_color, 
                     ribotype_027_color,  
                     ribotype_053_color, 
                     ribotype_078_color, 
                     ribotype_other_color)
  
  legend_names <- c("TreR L172I",
                    "TreR C171S", 
                    "Four Gene Insertion", 
                    "Severe Infection Outcome", 
                    "RT014", 
                    "RT015", 
                    "RT017", 
                    "RT027",
                    "RT053",
                    "RT078",
                    "Other Ribotype")
  legend_shapes <- c(empty_circle,
                     empty_circle, 
                     empty_circle, 
                     triangle,
                     circle, 
                     circle, 
                     circle, 
                     circle, 
                     circle, 
                     circle, 
                     circle)
  
  tree_cex <-  0.25

  # PLOT ----------------------------------------------------------------------#
  pdf(paper = "USr", 
      file = paste0("../figures/",
                    Sys.Date(), 
                    "_trehalose_tree.pdf"), 
      width = 10, 
      height = 8)
  plot(tree, 
       use.edge.length = TRUE, 
       label.offset = 0.0003, 
       type = "fan", 
       show.tip.label = show_labels, 
       cex = tree_cex, 
       no.margin = TRUE, 
       main = "", 
       rotate = 2.75)
  
  # ADD SCALE BAR
  add.scale.bar(font = 3)
  
  # ADD RIBOTYPE
  tiplabels(tip = c(1:Ntip(tree)),
            pch = circle, 
            cex = dot_size, 
            col = ribotype_other_color, 
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_027],
            pch = circle, 
            cex = dot_size, 
            col = ribotype_027_color,
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_014],
            pch = circle, 
            cex = dot_size, 
            col = ribotype_014_color, 
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_015],
            pch = circle, 
            cex = dot_size, 
            col = ribotype_015_color, 
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_017],
            pch = circle, 
            cex = dot_size,
            col = ribotype_017_color,
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_053],
            pch = circle, 
            cex = dot_size, 
            col = ribotype_053_color, 
            offset = ribotype_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% ribotype_078],
            pch = circle, 
            cex = dot_size, 
            col = ribotype_078_color, 
            offset = ribotype_offset)
  
  # ADD SEVERE OUTCOME
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% cases], 
            pch = triangle,
            cex = dot_size, 
            col = cases_color, 
            offset = cases_offset)
  
  # ADD TREHALOSE VARIANTS
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% four_gene_insertion], 
            pch = empty_circle, 
            cex = dot_size, 
            col = four_color, 
            offset = four_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% L172I], 
            pch = empty_circle,
            cex = dot_size, 
            col = L172I_color, 
            offset = L172I_offset)
  
  tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% C171S], 
            pch = empty_circle, 
            cex = dot_size, 
            col = C171S_color, 
            offset = C171S_offset)
  
  # ADD LEGEND
  legend("topright", 
         legend_names, 
         col = legend_colors, 
         pch = legend_shapes, 
         bty = "n", 
         title = "")
  dev.off()
} # end plot_trehalose_tree()