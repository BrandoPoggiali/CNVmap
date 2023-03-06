#' Plot single chromosome log2 coverages
#'
#' This function create a plot displaying log2 values in a specific chromosome. You can give 
#' to this function one or more genes in form of a vector and the function will plot them
#' base on their position and width in plot.This function#'can also plot a chromosome 
#' icon which is colored base on log2 values. It can have a color gradient which goes from red,
#' when there is a negative log2 value,to green when there is a positive log2 value.   
#'
#' @param cnr_data Dataframe containing .cnr file originated from the tool CNVkit.
#' @param cns_data Dataframe containing .cns file originated from the tool CNVkit.
#' @param chr Chromosome the function should plot The accepted codes are:
#' "chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
#' "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
#' "chr22", "chrX", "chrY". Default: "chr1"
#' @param chr_picture  Plot a cromosome icon colored base on log2 values. Default: FALSE
#' @param genes  A string with a gene or a vector containing gene names (HGNC ID). Default: NULL 
#' @param log2_line_col Color for the average log2 line. Default: "deepskyblue"
#' @param log2_threshold Numeric value to create a threshold line in the plot. Defaul: NULL 
#' @param log2_threshold_color Color for the log_threshold line. Default: "red"
#' @return A ggplot graph. 
#' @export
plot_single_chr <- function(cnr_data, cns_data, chr = "chr9", chr_picture = FALSE, genes = NULL, log2_line_col = "deepskyblue",
                            log2_threshold = NULL, log2_threshold_color = "red") {
  
  check_cnvkit_data(cnr_data, cns_data)
  
  if (!chr %in% available_chromosomes){
    stop(' You must supply one of these chromosomes code: "chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY". ', call. = FALSE) 
  }
  
  if (log2_line_col != "red" & length(log2_line_col) == 1) {
    if (!all(.areColors(log2_line_col))) {
      stop(paste("log2_line_col contains wrong color code:", log2_line_col), call. = FALSE)
    } 
  }
  
  if (!is.null(log2_threshold)){
    if (!is.numeric(log2_threshold)) {
      stop("log2_threshold parameter is not a number.", call. = FALSE)
    }
  }
  
  if (log2_threshold_color != "red" & length(log2_threshold_color) == 1) {
    if (!all(.areColors(log2_line_col))) {
      stop(paste("log2_line_col parameter contains wrong color code:", log2_threshold_color), call. = FALSE)
    } 
  }
  
  cnr_chr <- cnr_data[cnr_data$chromosome == chr,]
  cns_chr <- cns_data[cns_data$chromosome == chr,]
  plot_point <- (cnr_chr$end[1] - cnr_chr$start[1]) / 2
  centromere_position_chr <- centromere_positions[centromere_positions$CHROM == chr,]

    ###Plots------------------------------------------------------------------------------------------------------------

  graph <- ggplot() +
    geom_point(data = cnr_chr, aes(x = (start + plot_point) / 1000000, y = log2, pch = "."), color = "#616a6b", size = 0.8) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.minor.y = element_line(linewidth = 0.1, colour = "grey50", linetype = "dashed"),
          panel.grid.major.y = element_line(linewidth = 0.1, colour = "grey50", linetype = "dashed")) +
    #ggtitle(paste("Chromosome", stringr::str_remove(chr, "chr"))) + 
    xlab("Position (Mb)") + ylim(-3, 3) +
    geom_hline(yintercept = 0, linewidth=1.08 ) +  #set black line at y=0
    scale_x_continuous(expand = c(0, 0)) +  #set side layout graph
    geom_rect(data = as.data.frame(centromere_position_chr), inherit.aes = FALSE, aes(xmin = START / 1000000, xmax = END / 1000000, ymin = -Inf, ymax = Inf), color = "transparent", fill = "orange", alpha = 0.27) +  #Plot centromere position
    geom_segment(data = cns_chr, aes(x = start / 1000000, xend = end / 1000000, y = log2, yend = log2), lineend = "round", color = log2_line_col, size = 1.45) +
    theme(plot.title = element_text(hjust = 0.5, size = 22), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11), axis.title = element_text(size = 15), legend.position = "none")

  options(warn = -1)
  if (!is.null(log2_threshold)) {
    graph <- graph + geom_hline(yintercept = log2_threshold, linewidth = 0.35, color = "red", linetype="dashed") +
      scale_y_continuous(limits = c(-5, 2), breaks = c(2, 0, -2, -4, log2_threshold), minor_breaks = c(1, -1, -3, -5)) 
  }
  options(warn = 0)
  
  if (!is.null(genes)) {
    target_gene_info <- subset(all_genes[seqnames(all_genes) == gsub("chr", "", chr)], symbol %in% genes)
    
    if (length(target_gene_info) == 0) {
      warning(paste("No supplied gene/s were found in the chromosome", gsub("chr", "", chr)), call. = FALSE)
      } else {
    a <- seq(1.4, 3, by = 0.20)
    b <- seq(3.2, 5, by = 0.20)
    x <- b
    value_to_change <- "b"
    row=7
    for (row in 1:length(target_gene_info)) {
      gene_start <- start(target_gene_info[row,]) / 1000000
      gene_end <- end(target_gene_info[row,]) / 1000000
      gene_name <- mcols(target_gene_info[row, "symbol"])[1,1]
      
      #if vectors are empty fill in again
      if(length(x) == 1 && value_to_change == "a"){
        a <- seq(1.4,3, by = 0.20)
        a <- a[!a %in% last_value_a]
        x <- a
      }

      if(length(x)==1 && value_to_change=="b"){graph
        b <- seq(3.2,5,by=0.20)
        b <- b[! b %in% last_value_b]
        x <- b
      }

      y_value=sample(x,1)

      #remove value from vectors
      if (value_to_change=="a"){
        a=a[! a %in% y_value]
        last_value_a <- y_value
      }
      if (value_to_change=="b"){
        b=b[! b %in% y_value]
        last_value_b <- y_value
      }

      #change variable to choose the vector to use
      if (value_to_change == "a") {
        value_to_change <- "b"
        x <- b
      } else {
        value_to_change <- "a"
        x <- a
      }

      if (gene_end - gene_start <= 0.150){
      graph <- graph + geom_vline(xintercept = (gene_start + gene_end) / 2 , color="#f9a4a4", size=0.25 + gene_end - gene_start) +
               annotate("text", x = mean(gene_start, gene_end), y=-y_value, label = gene_name, col="red", size=2.80) +
        ylim(-5, 2)
      } else {
       graph <- graph + geom_rect(aes(xmin = !!gene_start, xmax = !!gene_end, ymin = -Inf , ymax = Inf), fill="#f9a4a4") +
         annotate("text", x = mean(gene_start, gene_end), y = -y_value,label = gene_name, col="red", size=2.80) + theme(legend.position = "none") +
         ylim(-5, 2)
    }}}}
    

#Plot
if (chr_picture){
  chromosome_length <- cnr_chr[nrow(cnr_chr),"end"] - cnr_chr[1,"start"]   
  
  centromere_space <- chromosome_length / 100 #Create a centromere which has a wisth of 1% of the total lenght of the chromosome
  centromere <- rowMeans(centromere_position_chr[,c(2,3)])
  
  color_cns_chr <- cns_chr
  color_cns_chr$log2 <- ifelse(color_cns_chr$log2 < -2, -2, color_cns_chr$log2)
  color_cns_chr$log2 <- ifelse(color_cns_chr$log2 > 2, 2, color_cns_chr$log2)
  
  start_centr <- centromere - centromere_space
  end_centr <- centromere + centromere_space
  
  color_cns_chr[nrow(color_cns_chr) + 1,] <- c(chr, start_centr, end_centr, NA, NA, NA, NA, NA, NA, NA)
  color_cns_chr$start <- as.numeric(color_cns_chr$start) 
  color_cns_chr$end <- as.numeric(color_cns_chr$end) 
  color_cns_chr$log2 <- as.numeric(color_cns_chr$log2)
  
  
  chr_rect <- ggplot(color_cns_chr) +
    geom_rect(aes(xmin = start / 1000000, xmax = end / 1000000, ymin = 1, ymax = 2, fill = log2)) +
    scale_x_continuous(expand = c(0, 0)) + scale_fill_gradient2(low = "#d20e07", mid = "#bababa", high = "#18852f",
                                                                        midpoint = 0, na.value = "white") +
    theme_void() + theme(legend.position = "none") + ggtitle(paste("Chromosome", stringr::str_remove(chr, "chr"))) +
    theme(plot.title = element_text(hjust = 0.5, size = 22))
  
  cent_df <- data.frame(x=c(start_centr / 1000000, #build centromere dataframe for plot
                         start_centr / 1000000,
                         centromere / 1000000,
                         end_centr / 1000000,
                         end_centr / 1000000,
                         centromere / 1000000),
                     y=c(1,2,1.5,  1,2,1.5))
  
  #cent_df <- .create_centromere(centromere, centromere_space)
  chr_rect <- chr_rect + geom_polygon(data=cent_df, mapping = aes(x = x, y = y), fill="black") 
  
  #Fill empty spaces
  centromere_loc <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chrom = chr, start = cent_df[1,1] * 1000000, end = cent_df[5,1] * 1000000))
  empty <- GenomicRanges::makeGRangesFromDataFrame(color_cns_chr)
  empty <- append(empty, centromere_loc)
  empty <- as.data.frame(GenomicRanges::reduce(empty))
  empty <- data.frame(start = empty$end[-length(empty$end)], end = empty$start[-1])
  empty$start <-  empty$start / 1000000
  empty$end <-  empty$end / 1000000    
  
  
  
  chr_rect <- chr_rect + geom_segment(data = empty, aes(x = start, y = 1.01, xend = end, yend = 1.01), color = "#bababa") +
    geom_segment(data = empty, aes(x = start, y = 1.99, xend = end, yend = 1.99), color = "#bababa")
    #geom_segment(data = empty, x = start, y = 2, xend = end, yend = 2)
  
  
  return(ggpubr::ggarrange(chr_rect, graph, heights = c(0.14, 0.86), ncol=1, nrow=2, align="v"))
} else {
  graph <- graph + ggtitle(paste("Chromosome", stringr::str_remove(chr, "chr")))
  return(graph)
}}


.create_centromere <- function(centromere_location, dim_per_side){
  start_centr <- (centromere_location - dim_per_side) / 1000000
  end_centr <- (centromere_location + dim_per_side) / 1000000
  #abs(start_centr - end_centr)/2
  dim_rect <- ((dim_per_side * 2) / 1000000) / 5000 #5,000 is the value to obtain the centromere structure
  
  cent <- data.frame(start = seq(start_centr + dim_rect, end_centr, dim_rect),
                     end=seq(start_centr + (dim_rect*2), end_centr + dim_rect, dim_rect),
                     min=seq(1.00005,2,0.0002),
                     max=seq(2,1.00005,-0.0002))
  return(cent)
}


