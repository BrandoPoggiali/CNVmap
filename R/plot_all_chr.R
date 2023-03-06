#' Plot all chromosome log2 coverages
#'
#' This function create a plot displaying log2 ratios in all chromosomes. 
#'
#' @param cnr_data Dataframe containing .cnr file originated from the tool CNVkit.
#' @param cns_data Dataframe containing .cns file originated from the tool CNVkit.
#' @param only_autosomal Parameter to plot only autosomal chromosomes. Default: FALSE
#' @param chr_text_size  Size of the chromosome number text displayed at the top of the plot. Default: 4
#' If you don't want any text set chr_text_size = NULL
#' @param chr_colors A vector of colors. Colors will be used om the plot  base on their order in the vector. I you plot only autosomes you must input 22 colors,
#' if you plot also sex chromosomes you must input 24 colors. Default: NULL
#' @param log2_line_col Color for the average log2 line. Default: "red"
#' @return A ggplot graph.
#' @export
plot_all_chr <- function(cnr_data, cns_data, only_autosomal = FALSE, chr_text_size = 4, chr_colors = NULL, log2_line_col = "red"){
  
  check_cnvkit_data(cnr_data, cns_data)
                    
  if (!is.numeric(chr_text_size)  & length(chr_text_size) != 1){
    if (is.null(chr_text_size)) {
      chr_text_size <- 0 } else {
    stop(paste(substitute(chr_text_size),"should be a one numeric value!"),call.=FALSE)
  }}
   
  if (!is.logical(only_autosomal)) {
    stop(paste(only_autosomal," must be a logical value!"),call.=FALSE)
  }
  
  if (!is.null(chr_colors)) {
    if (only_autosomal == TRUE & length(chr_colors) < 22) stop(paste(substitute(chr_colors),"contains less than 22 colors"), call. = FALSE)
    if (only_autosomal == FALSE & length(chr_colors) < 24) stop(paste(substitute(chr_colors),"contains less than 24 colors"), call. = FALSE)
    if (!all(.areColors(chr_colors))) {
      not_colors <- which(areColors(chr_colors) %in% FALSE)
      stop(paste(substitute(chr_colors),"contains non color code/s, check element/s:", toString(not_colors)), call. = FALSE)
    }} else {
      #Default colors
      chr_colors <- c("#FDB462","#7FC97F","#CAB2D6","#84adff","#FFF2AE","#E6AB02","#A6D854","#d4dd1d","#FB9A99","#999999","#ac5bb5","#ceb047",
                      "#BF5B17","#1B9E77","#e69fe2","#FB8072","#72e3f1","#A6761D","#FDCDAC","#ffd680","#E78AC3","#2365a3","#2a9239", "#ff9d79")
    }
    
  if (log2_line_col != "red" & length(log2_line_col) == 1) {
    if (!all(.areColors(log2_line_col))) {
      stop(paste(substitute(log2_line_col),"contains wrong color code:", log2_line_col), call. = FALSE)
    } 
  }
  
  
  pb <- txtProgressBar(min = 0, max = 28, style = 3, width = 50, char = "=") #create bar loading
  cnr_data$start <- as.numeric(cnr_data$start) #give correct structure for the analysis
  cnr_data$end <- as.numeric(cnr_data$end)
  cns_data$start <- as.numeric(cns_data$start)
  cns_data$end <- as.numeric(cns_data$end)
  plot_point <- (cnr_data$end[1]-cnr_data$start[1])/2 
  
  setTxtProgressBar(pb, 1)
  ###Elaborate dataframes for plotting-------------------------------------------------------------------------------------------------
  all_chr <- paste0("chr",c(1:22,"X","Y"))

  #Save chromosomes lengths
  chromosome_lengths <- c()
  for (c in all_chr){
    chromosome_lengths <- append(chromosome_lengths,tail(cnr_data[cnr_data$chromosome == c ,],n=1)[1,3])
  }
  
  #Elaborating cnr_data for plot
  n_loading_bar <- 1
  
  for (c in all_chr){
  
    if (c != "chrY"){
    next_chr <- all_chr[which(all_chr == c) + 1]
  
    value=tail(cnr_data[cnr_data$chromosome == c ,],n=1)[1,3] #take last value chromosome position
  
    cnr_data$start <- ifelse(cnr_data$chromosome == next_chr,cnr_data$start + value,cnr_data$start)
    cnr_data$end <- ifelse(cnr_data$chromosome == next_chr, cnr_data$end + value,cnr_data$end)
    setTxtProgressBar(pb, n_loading_bar)
    n_loading_bar <- n_loading_bar + 0.5
    }}
  
  
  #Elaborating cns_data data for plot
  for (c in all_chr){
  
    if (c != "chrY"){
      next_chr <- all_chr[which(all_chr == c) + 1]
      value=tail(cnr_data[cnr_data$chromosome == c ,],n=1)[1,3]
  
      cns_data$start <- ifelse(cns_data$chromosome == next_chr, cns_data$start + value,cns_data$start)
      cns_data$end <- ifelse(cns_data$chromosome == next_chr, cns_data$end + value,cns_data$end)
      setTxtProgressBar(pb, n_loading_bar)
      n_loading_bar <- n_loading_bar + 0.5
    }}
  
  #Creating dataframe with separating line
  sep_chrom <- data.frame(matrix(ncol=1,nrow=24))
  colnames(sep_chrom) <- "position"
  n=1
  for (c in all_chr){
  
      sep_chrom[n,1] <- tail(cnr_data[cnr_data$chromosome == c ,],n=1)[1,3]
      n=n+1
    }
  setTxtProgressBar(pb, n_loading_bar+0.5)
  sep_chrom <- cbind(chromosome=all_chr,sep_chrom)
  auto_chromosome_lenghts <- as.numeric(chromosome_lengths)
  sep_chrom$position_text <- sep_chrom$position - (auto_chromosome_lenghts/2)
  sep_chrom$position_text[nrow(sep_chrom)-1] <- sep_chrom$position_text[nrow(sep_chrom)-1] + 1300000
  sep_chrom$position_text[nrow(sep_chrom)] <- sep_chrom$position_text[nrow(sep_chrom)] + 4000000
  setTxtProgressBar(pb, n_loading_bar+1)
  ###Plots------------------------------------------------------------------------------------------------------------
  
  names(chr_colors) <- all_chr
  if (only_autosomal){
    cnr_data <- cnr_data[!cnr_data$chromosome == "chrX" & !cnr_data$chromosome == "chrY",]
    cns_data <- cns_data[!cns_data$chromosome == "chrX" & !cns_data$chromosome == "chrY",]
    chr_colors <- head(chr_colors,-2) 
    names(chr_colors) <- head(all_chr,-2)
    sep_chrom <- head(sep_chrom,-2)
  } 
  setTxtProgressBar(pb, n_loading_bar+2)
    graph <- ggplot() +
    geom_point(data=cnr_data, aes(x=(start+plot_point)/1000000, y=log2, color=chromosome, pch="."), size = 0.8) +
    theme(legend.position="none", panel.background = element_rect(fill = "white", colour = "#d2d3d3"), #panel.grid.minor.y = element_line(size=0.1, colour="#d2d3d3", linetype="dashed"),
          panel.grid.major.y = element_line(size=0.1, colour = "grey50",linetype="dashed")) +
    scale_color_manual(values = chr_colors) +
    ylim(-2.5,2.5) + xlab("Chromosomes") +
    geom_hline(yintercept = 0, size=1.08 ) +  #set black line at y=0
    scale_x_continuous( expand = c(0.002, 0.0008)) +
      # geom_rect(data=as.data.frame(centromere_positions), inherit.aes=FALSE, aes(xmin=START/1000000, xmax=END/1000000, ymin=-Inf,ymax=Inf ), color="transparent", fill="orange", alpha=0.27) + #Plot centromere position
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    geom_segment(data=cns_data, aes(x=start/1000000, xend=end/1000000, y=log2, yend=log2), lineend="round", color= log2_line_col, size=1.35) +
    geom_vline(xintercept = head(sep_chrom$position,-1)/1000000) +
    annotate("text",x=sep_chrom$position_text/1000000, y=2.2 , label= gsub("chr","",sep_chrom$chromosome), size=chr_text_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_text(face="bold",size=15) ,axis.title.y = element_text(face="bold",size=12.5), axis.text.y = element_text(size = 11))
  
  setTxtProgressBar(pb, 27)
  options(warn=-1)
  print(graph)
  options(warn=0)
  Sys.sleep(10)
  setTxtProgressBar(pb, 28)
  close(pb)
  
  }


#Test if color codes are valid
.areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

