#' save_four_gene_insertion
#' # 1. Read in the pan genome file. 
# 2. Extract the sequences for the four gene trehalose insertion. 
# 3. Save in a new fasta file. 
#' @param pan_genome_path 
#'
#' @return
#' @export
#'
#' @examples
save_four_gene_insertion <- function(pan_genome_path){
  # Note: command: roary -p 8 -r -e -n -v -i 70 *.gff
  
  # Based on BLAST results these are the four genes of interest:  
  # Name  | Roary name  
  # treA2 | treA_1  
  # ptsT  | treB  
  # treX  | group_5209  
  # treR2 | treR_1  
  
  # 1. Read in the pan genome file. 
  pan_genome_seq <- 
    seqinr::read.fasta(file = pan_genome_path, 
                       seqtype = "DNA", 
                       as.string = TRUE, 
                       forceDNAtolower = FALSE)
  
  # 2. Extract the sequences for the four gene trehalose insertion. 
  trehalose_indices <- rep(0, 4)
  trehalose_indices[1] <- grep("treA_1$", seqinr::getAnnot(pan_genome_seq))
  trehalose_indices[2] <- grep("treB$", seqinr::getAnnot(pan_genome_seq))
  trehalose_indices[3] <- grep("group_5209$", seqinr::getAnnot(pan_genome_seq))
  trehalose_indices[4] <- grep("treR_1$", seqinr::getAnnot(pan_genome_seq))
  
  four_gene_trehalose_insertion <- pan_genome_seq[trehalose_indices]
  seq_names <- c("treA2", "ptsT", "treX", "treR2")
  
  # 3. Save in a new fasta file.  
  seqinr::write.fasta(sequences = four_gene_trehalose_insertion, 
                      names = seq_names, 
                      file.out = 
                        "../data/outputs/four_gene_trehalose_insertion.fna", 
                      open = "w", 
                      as.string = TRUE, 
                      nbchar = 60)
  
  seqinr::write.fasta(sequences = four_gene_trehalose_insertion, 
                      names = seq_names, 
                      file.out = 
                        "../data/outputs/supplementary_file_1.fna", 
                      open = "w", 
                      as.string = TRUE, 
                      nbchar = 60)
} # end save_four_gene_insertion()