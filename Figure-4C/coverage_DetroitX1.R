## Composite plots for DetroxitX1 cell lines
## Create composite plots based on
library(smoother)
library(bamsignals)
library(GenomicRanges)
library(ggplot2)

composite_plot <- function(cdtdir, window = 0.1, normalize = TRUE){
  #composite_plot <- function(cdtdir, window = 0.1, normalize = TRUE, factors){
  
  df = data.frame(distance=integer(),count=double(),samples=character())
  for(cdt in list.files(cdtdir)){
    if(tools::file_ext(cdt) == "cdt"){
      
      cdtmat = as.matrix(read.table(file.path(cdtdir,cdt),header=TRUE, sep = "\t", as.is = TRUE, row.names = 1, check.names = FALSE))
      #cdtmat = as.matrix(read.table(file.path(cdtdir,cdt),header=TRUE, sep = "\t", as.is = TRUE, check.names = FALSE))
      length = dim(cdtmat)[2]
      summation = colSums(cdtmat)
      
      cdtmat = apply(cdtmat, c(1,2), function(x) {(as.integer(x))})
      #summation = colSums(cdtmat[,-which(colnames(cdtmat) == "GENE")])
      ##summation = summation*factors[cdt,]
      #group = tools::file_path_sans_ext(basename(cdt))
      group = unlist(strsplit(tools::file_path_sans_ext(basename(cdt)),"_"))[1]
      y = smth(summation,window = window,method = "gaussian")
      # Divide by the max to plot all curves on the same graph
      if(normalize){
        y = y/max(na.omit(y))
        
      }
      
      tmp = data.frame(distance = as.integer(names(summation)), count = y, samples = group)
      df = rbind(df,tmp)
      
      
    }
  }
  if(normalize){
    ylim = 1
  } else{
    ylim = max(na.omit(df$count))
  }
  ggplot(df, aes(x=distance, y=count,group=samples, color=samples)) + geom_line(size=0.75) + lims(y=c(0,ylim)) + labs(y="Average Counts",x="Distance from ref (bp)") + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                                                                                                          panel.grid.major = element_blank(),
                                                                                                                                                                                          panel.grid.minor = element_blank(),
                                                                                                                                                                                          panel.background = element_blank())
}
## Plot the agregated reads along TEAD DMSO peaks
pdf("DetroitX1_TEAD.pdf")
composite_plot('bam/norm_CDT_NGS2957/', window = 0.1, normalize = FALSE)
dev.off()

