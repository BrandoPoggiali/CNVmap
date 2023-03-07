#' Plot gene log2 coverage
#'
#' 
#' This function create a plot displaying log2 values in a specific chromosome. You can give 
#' to this function one or more genes in form of a vector and the function will plot them
#' base on their position and width in plot.This function#'can also plot a chromosome 
#' icon which is colored base on log2 values. It can have a color gradient which goes from red,
#' when there is a negative log2 value,to green when there is a positive log2 value.   
#'
#' @param cnr_data Dataframe containing .cnr file originated from the tool CNVkit.
#' @param cns_data Dataframe containing .cns file originated from the tool CNVkit.
#' @param genes  A string with a gene name (HGNC ID). Default: "TP53"
#' @param all.transcripts Boolean value for make the plot showing all transcripts of the target gene. Default = FALSE 
#' @param gene_structure Boolean value for make the plot showing canonical transcript of the target gene. Default = TRUE 
#' @param regulatory.elements Boolean value for make the plot displaying regulatory elements from ... Default: FALSE
#' @param log2_line_col Color for the average log2 line. Default: "deepskyblue"
#' @param log2_threshold Numeric value to create a threshold line in the plot. Defaul: NULL
#' @param log2_threshold_color Color for the log_threshold line. Default: "red"
#' @return A matrix of the infile
#' @export
gene = "ARHGEF10"
gene = "PTPRD"
plot_single_gene <- function(cnr_data, cns_data, gene = "TP53", all.transcripts = FALSE, gene_structure = TRUE,
                             regulatory.elements = FALSE, log2_line_col = "deepskyblue", log2_threshold = NULL,
                             log2_threshold_color = "red"){
  
  check_cnvkit_data(cnr_data, cns_data)
  
  if (log2_line_col != "red" & length(log2_line_col) == 1) {
    if (!all(.areColors(log2_line_col))) {
      stop(paste("log2_line_col parameter contains wrong color code:", log2_line_col), call. = FALSE)
    } 
  }
  
  if (!is.logical(all.transcripts)) stop("all.transcripts must be a boolean value", call. = FALSE)
  if (!is.logical(gene_structure)) stop("gene_structure must be a boolean value", call. = FALSE)
  
  cnn_trs_id <- all_genes[all_genes$symbol == gene ]$canonical_transcript
  if (length(cnn_trs_id) == 0) {
    stop(paste(substitute(gene),"is not a valid gene name."), call. = FALSE)
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
  
  cnn_trs_id <- all_genes[all_genes$symbol == gene][1]$canonical_transcript
  trs_start <- start(all_transcripts[names(all_transcripts) == cnn_trs_id])
  trs_end <- end(all_transcripts[names(all_transcripts) == cnn_trs_id])
  gene_chr <- as.vector(seqnames(all_genes[all_genes$symbol == gene][1]))
  gene_strand <- as.vector(strand(all_genes[all_genes$symbol == gene])[1])
  gene_exons <- as.data.frame(all_exons[all_exons$tx_id == cnn_trs_id])
  
  #Filter dataset (I should transform cns and cnr file in Granges object and make the filtering using Granges FindOverlap)
  cnr_data <- makeGRangesFromDataFrame(cnr_data, keep.extra.columns = TRUE)
  cns_data <- makeGRangesFromDataFrame(cns_data, keep.extra.columns = TRUE)
  trs_cnn_gr <- all_transcripts[names(all_transcripts) == cnn_trs_id]

  x <- seqlevels(trs_cnn_gr)
  newStyle <- mapSeqlevels(x,"UCSC")
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  trs_cnn_gr <- GenomeInfoDb::renameSeqlevels(trs_cnn_gr, newStyle)
  
  cnr_gene <- as.data.frame(subsetByOverlaps(cnr_data, trs_cnn_gr))
  cns_gene <- as.data.frame(subsetByOverlaps(cns_data, trs_cnn_gr))

  #Code to display the transcript with exons and introns. 
  if (!all.transcripts){
    arrow_tick_trs <- seq(trs_start, trs_end, by = abs(trs_start - trs_end) / 40)
    
    if (gene_strand == "+"){
      arrow_start <- head(arrow_tick_trs, -1)
      arrow_end <- arrow_tick_trs[-1]
      no_arrows_areas <- data.frame(start = gene_exons$start ,
                                    end = gene_exons$end + (1 * (trs_end - trs_start) / 100))
    } else {
      arrow_start <- arrow_tick_trs[-1]
      arrow_end <- head(arrow_tick_trs, -1)
      no_arrows_areas <- data.frame(start = gene_exons$start - (1 * (trs_end - trs_start) / 100),
                                    end = gene_exons$end)
    }  
    
    graph_start <- ifelse(cnr_gene$start[1] < GenomicRanges::start(trs_cnn_gr), cnr_gene$start[1], GenomicRanges::start(trs_cnn_gr))
    graph_end <- ifelse(cnr_gene$end[nrow(cnr_gene)] > GenomicRanges::end(trs_cnn_gr), cnr_gene$end[nrow(cnr_gene)], GenomicRanges::end(trs_cnn_gr))
  
    graph <- ggplot() +
      theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.minor.y = element_line(linewidth = 0.2, colour = "#aaaaaa", linetype = "dashed"),
          panel.grid.major.y = element_line(linewidth = 0.1, colour = "grey50", linetype = "dashed"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
      geom_segment(data = cnr_gene, aes(x = start / 1000000, xend = end / 1000000 , y = log2, yend = log2), lineend = "round", color = "#616a6b", size=1.5) +
      ggtitle(paste0(cnr_gene$seqnames[1], ": ", gene, " gene")) + xlab("Position (Mb)") + #ylim(-4,2) +
      theme(plot.title = element_text(hjust = 0.5, size = 22)) +
      geom_hline(yintercept = 0, linewidth = 1.08 ) +  #set black line at y=0
      coord_cartesian(xlim = c(graph_start / 1000000, graph_end / 1000000)) +
      geom_segment(aes(x = cns_gene$start / 1000000, xend = cns_gene$end / 1000000, y=cns_gene$log2, yend=cns_gene$log2), lineend = "round"  ,color = log2_line_col, linewidth = 1.25) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11), axis.title = element_text(size = 14)) 
  
  if (!gene_structure) {
    graph <- graph + scale_y_continuous(limits = c(-2.5,2.5), expand = c(0.01, 0.01))
  } else {
    
    #Aesthetics for wide and narrow exons
    wide_gene_exons <- gene_exons[((gene_exons$end - gene_exons$start) * 100) / (trs_end - trs_start) > 0.2,]
    narrow_gene_exons <- gene_exons[!((gene_exons$end - gene_exons$start) * 100) / (trs_end - trs_start) > 0.2,]
    narrow_gene_exons$exon_center <- rowMeans(narrow_gene_exons[, c("start", "end")])  
    
    #Aesthetics for overlapping arrows with exons
    
    index_to_remove <- c()
    n=1
    for (arrowhead in arrow_end) {
      overlap <- nrow(no_arrows_areas[no_arrows_areas$end >= arrowhead & no_arrows_areas$start <= arrowhead,])
      if (overlap > 0) {
        index_to_remove <- c(index_to_remove, n)
      }
      n=n+1
    }
    
    arrow_start <- arrow_start[-index_to_remove]
    arrow_end <-  arrow_end[-index_to_remove]
    
    graph <- graph + geom_segment(aes(x = arrow_start / 1000000, xend = arrow_end / 1000000, y = -3.5, yend = -3.5), color = "#3b7fc6", linewidth = 0.5, arrow = arrow(length = unit(0.1, "inches"))) +
      geom_segment(aes(x = trs_start / 1000000, xend = trs_end / 1000000, y = -3.5, yend = -3.5), color = "#3b7fc6", linewidth = 0.5) +
      geom_rect(data = wide_gene_exons, aes(xmin = start / 1000000, xmax =end/ 1000000, ymin = -3.7, ymax = -3.30), fill = "#1c2898") +
      geom_segment(data = narrow_gene_exons, aes(x = exon_center / 1000000, xend = exon_center / 1000000, y = -3.7, yend = -3.3), color = "#1c2898") +
      scale_y_continuous(limits = c(-4,2.5), expand = c(0.01, 0.01))
   
  }} else {
    gene_id <- names(all_genes[all_genes$symbol == gene ])
    all_exons_gene <- all_exons[all_exons$gene_id == gene_id]
    gene_transcripts <- unique(all_exons_gene$tx_id)
    all_transcripts_gene <- all_transcripts[all_transcripts$gene_id == gene_id]
    n_all_transcripts <- length(gene_transcripts)
    
    #Subset data
    start_gene_region <- min(start(all_transcripts_gene))
    end_gene_region <- max(end(all_transcripts_gene))
    
    cnr_data_trascr <- plyranges::filter(cnr_data, seqnames == paste0("chr", gene_chr) & end >=  start_gene_region & start <=  end_gene_region)
    cnr_data_trascr <- as.data.frame(cnr_data_trascr)
    cns_data_trascr <- plyranges::filter(cns_data, seqnames == paste0("chr", gene_chr) & end >=  start_gene_region & start <=  end_gene_region)
    cns_data_trascr <- as.data.frame(cns_data_trascr)
    
    
    if (regulatory.elements & n_all_transcripts > 4){
       y_value <- -4.6
       y_value_line <- -4.6
    } else {
       y_value <- -4.5
       y_value_line <- -4.5 
    }
    
    graph <- ggplot() +
    theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.1, colour = "grey50", linetype = "dashed"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_segment(data = cnr_data_trascr, aes(x = start / 1000000, xend = end / 1000000 , y = log2, yend = log2), lineend = "round", color = "#616a6b", size=1.7) +
    ggtitle(paste0(cnr_data_trascr$seqnames[1], ": ", gene, " gene")) + xlab("Position (Mb)") + #ylim(-4,2) +
    geom_hline(yintercept = 0, linewidth = 1.08 ) +  #set black line at y=0
    coord_cartesian(xlim = c(cnr_data_trascr$start[1] / 1000000, cnr_data_trascr$end[nrow(cnr_data_trascr)] / 1000000)) +
    geom_segment(data = cns_data_trascr, aes(x = start / 1000000, xend = end / 1000000, y = log2, yend = log2), lineend = "round", color = log2_line_col, linewidth = 1.4) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 22)) +
    scale_y_continuous(limits = c(y_value, 2), expand = c(0.01, 0.01), breaks = c(2, 1, 0, -1, -2)) +
    scale_x_continuous(expand = c(0, 0))
    
    
    arrow_tick <- seq(start_gene_region, end_gene_region, by = abs(start_gene_region - end_gene_region) / 40)
    
    #calculate value for plotting all the transcripts
    
    if (n_all_transcripts > 40){
      warning(paste(gene, "has", n_all_transcripts, "transcripts; only 40 trascripts are plotted for graphic reason!"), call. = FALSE)
      gene_transcripts <- gene_transcripts[1:40] 
      n_all_transcripts <- 40
      }
    
    
    exon_height <- 0.2 #exon_height
    intron_height <- 0.5
    arrow_size <- 0.1
    #n_all_transcripts <- 20
    #Set value for transcript aesthetics
    if (n_all_transcripts > 4){
      aesthetics_parameter <- ((15 - 1) / (80 - 4)) * (n_all_transcripts - 4) + 1
      aesthetics_intron_parameter <- ((3.5 - 1) / (40 - 4)) * (n_all_transcripts - 4) + 1
      aesthetics_arrow_parameter <- ((2.5 - 1) / (40 - 4)) * (n_all_transcripts - 4) + 1
    } else {
      aesthetics_parameter <- 1
      aesthetics_intron_parameter <- 1
      aesthetics_arrow_parameter <- 1
      }
    
    exon_height <- exon_height / aesthetics_parameter
    intron_height <- (intron_height / (aesthetics_parameter / ((n_all_transcripts * aesthetics_intron_parameter) / n_all_transcripts))) 
    arrow_size <- (arrow_size / (aesthetics_parameter / ((n_all_transcripts * aesthetics_arrow_parameter) / n_all_transcripts))) 
    
    
    
    n <- 1
    #transcript <- gene_transcripts[n]
    for (transcript in gene_transcripts){
      trs_start <- start(all_transcripts_gene[all_transcripts_gene$tx_id == transcript])
      trs_end <- end(all_transcripts_gene[all_transcripts_gene$tx_id == transcript])
      all_exons_gene_trascr <- as.data.frame(all_exons_gene[all_exons_gene$tx_id == transcript])
      
      arrow_tick_trs <- arrow_tick[arrow_tick >= trs_start & arrow_tick <= trs_end]
      
      
      if (gene_strand == "+"){
        arrow_start <- head(arrow_tick_trs, -1)
        arrow_end <- arrow_tick_trs[-1]
        no_arrows_areas <- data.frame(start = all_exons_gene_trascr$start ,
                                      end = all_exons_gene_trascr$end + (1 * (trs_end - trs_start) / 100))
      } else {
        arrow_start <- arrow_tick_trs[-1]
        arrow_end <- head(arrow_tick_trs, -1)
        no_arrows_areas <- data.frame(start = all_exons_gene_trascr$start - (1 * (trs_end - trs_start) / 100),
                                      end = all_exons_gene_trascr$end)
      }  
      
      
      #Aesthetics for wide and narrow exons
      wide_gene_exons_trascr <- all_exons_gene_trascr[((all_exons_gene_trascr$end - all_exons_gene_trascr$start) * 100) / (end_gene_region - start_gene_region) > 0.2,]
      narrow_gene_exons_trascr <- all_exons_gene_trascr[!((all_exons_gene_trascr$end - all_exons_gene_trascr$start) * 100) / (end_gene_region - start_gene_region) > 0.2,]
      narrow_gene_exons_trascr$exon_center <- rowMeans(narrow_gene_exons_trascr[, c("start", "end")])  
      
      #Aesthetics for overlapping arrows with exons
      index_to_remove <- c()
      v=1
      for (arrowhead in arrow_end) {
        overlap <- nrow(no_arrows_areas[no_arrows_areas$end >= arrowhead & no_arrows_areas$start <= arrowhead,])
        if (overlap > 0) {
          index_to_remove <- c(index_to_remove, v)
        }
        v=v+1
      }
      
      if (length(arrow_tick_trs) != 0 & !is.null(index_to_remove)){
        arrow_start <- arrow_start[-index_to_remove]
        arrow_end <-  arrow_end[-index_to_remove]
      }
      #y_value <- position_transcripts[n]
      y_value <- y_value + (exon_height * 2) + (exon_height / 3)
      
      graph <- graph + geom_segment(aes(x = !!arrow_start / 1000000, xend = !!arrow_end / 1000000, y = !!y_value, yend = !!y_value), color = "#3b7fc6", linewidth = intron_height, arrow = arrow(length = unit(arrow_size, "inches"))) +
        geom_segment(aes(x = !!trs_start / 1000000, xend = !!trs_end / 1000000, y = !!y_value, yend = !!y_value), color = "#3b7fc6", linewidth = intron_height) +
        geom_rect(aes(xmin = !!wide_gene_exons_trascr$start / 1000000, xmax = !!wide_gene_exons_trascr$end/ 1000000, ymin = !!y_value - exon_height, ymax = !!y_value + exon_height), fill = "#1c2898") +
        geom_segment(aes(x = !!narrow_gene_exons_trascr$exon_center / 1000000, xend = !!narrow_gene_exons_trascr$exon_center / 1000000, y = !!y_value - exon_height, yend = !!y_value + exon_height), color = "#1c2898")
      n <- n + 1
    }}
  
  
   
  if (regulatory.elements & all.transcripts){
    reg_elements_chr <- regulatory_elements %>% plyranges::filter_by_overlaps(as(paste0(paste0("chr", gene_chr), ":", cnr_data_trascr$start[1], "-", cnr_data_trascr$end[nrow(cnr_data_trascr)]), 'GRanges'))
    reg_elements_chr <- as.data.frame(reg_elements_chr)
    #convert 0-based regulatory element data frame in 1-based
    reg_elements_chr$start <- reg_elements_chr$start +1
    y_value <- y_value + (exon_height * 2) + (exon_height / 3)
    
    graph <- graph + geom_rect(data=reg_elements_chr, aes(xmin = start / 1000000, xmax = end / 1000000, ymin = !!y_value - exon_height, ymax = !!y_value + exon_height, fill = ucscLabel)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("prom" = "red", "enhP" = "orange", "enhD" = "yellow", "K4m3" = "pink", "CTCF" = "blue")) 
  }
  
  if (regulatory.elements & !all.transcripts & gene_structure){
    reg_elements_chr <- regulatory_elements %>% plyranges::filter_by_overlaps(as(paste0(paste0("chr", gene_chr), ":", trs_start, "-", trs_end), 'GRanges'))
    reg_elements_chr <- as.data.frame(reg_elements_chr)
    #convert 0-based regulatory element data frame in 1-based
    reg_elements_chr$start <- reg_elements_chr$start +1
   
    graph <- graph + geom_rect(data=reg_elements_chr, aes(xmin = start / 1000000, xmax = end / 1000000, ymin = -3.1, ymax = -2.7, fill = ucscLabel)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("prom" = "red", "enhP" = "orange", "enhD" = "#ffe02c", "K4m3" = "pink", "CTCF" = "#40c8f8")) +
      scale_y_continuous(limits = c(-4, 2.5), breaks = c(2, 0, -2, -4, log2_threshold), minor_breaks = c(1, -1))
  }
  
  if (!is.null(log2_threshold)) {
    graph <- graph + geom_hline(yintercept = log2_threshold, linewidth = 0.25, color = "red", linetype="dashed") +
      scale_y_continuous(limits = c(y_value_line, 2.5), expand = c(0.01, 0.01), breaks = c(2, 1, 0, -1, -2, log2_threshold)) 
      
  }
  
  
  return(graph) 
}






