#' Generating Volcano plot for differential gene expression results
#' 
#' Function uses result gene list from differential expression analysis to 
#' prepare volcano plot.
#' 
#' Input gene list is filtered for differentially expressed genes according to given significance thresholds
#' for p-value and log foldchange. Significance thresholds are indicated as dashed lines in the plot.
#' Up and down regulated genes can be colored accordingly. The counts of significantly differentially expressed
#' genes is given in the upper plot margin.
#' 
#' @param data data.frame containing gene ID or Symbol (optional), p-value and logarithmized foldchange.
#' Column names can be given in \code{column_name_geme}, \code{column_name_p} and \code{column_name_fc}, respectively.
#' @param column_name_geme character with optional column name of gene names or symbols in \code{data}.
#' @param column_name_p character with column name of p-values in \code{data}.
#' @param column_name_fc character with column names of logarithmized fold changes in \code{data}.
#' @param title character with title of volcano plot.
#' @param p_value_threshold numeric significance threshold for p-value.
#' @param fc_threshold numeric significance threshold for logarithmized fold changes.
#' @param xlabel character with label for x-axis
#' @param ylabel character with label for y-axis
#' @param color_up_reg character with color name for indicating up regulated genes.
#' @param color_down_reg character with color name for indicating down regulated genes.
#' @param color_not_reg character with color name for indicating genes not differentially regulated.
#' 
#' @return dataframe with filtered gene list. Volcano plot is generated as side effect. 
#' 
#' @author Frank Ruehle
#' 
#' @export   



plot_volcano <- function (data, column_name_geme = NULL, column_name_p ="padj", column_name_fc = "log2FoldChange", 
                          title = "Volcano plot", p_value_threshold = 0.05, fc_threshold = log2(1.5),
                          xlabel="log fold change", ylabel="-log10 p-value (adjusted)", 
                          color_down_reg= "red", color_up_reg= "darkgreen", color_not_reg = "darkgray")
                          {


  # filter for significance
  data_filt  <- filterGeneLists(data,
                                newheader=NULL,
                                filtercat1 = column_name_p,
                                filtercat1.decreasing = FALSE,
                                filtercat1.function = identity,
                                filtercat1.threshold= p_value_threshold,
                                filtercat2 = column_name_fc,
                                filtercat2.decreasing = TRUE,
                                filtercat2.function = abs,
                                filtercat2.threshold = fc_threshold)
  
  data_filt$down <- data_filt[, column_name_fc] < 0
  data_filt$up <- data_filt[, column_name_fc] > 0

  plot(data[,column_name_fc], -log10(data[,column_name_p]), main = title,
     col= color_not_reg, pch=16, cex=0.7, xlab=xlabel, ylab=ylabel) 

  if(!is.null(p_value_threshold)) {abline(h=-log10(p_value_threshold),lty=2)}

  if(!is.null(fc_threshold)) {
    abline(v=c(-fc_threshold, fc_threshold), lty=2, col = c(color_down_reg, color_up_reg))

        if(nrow(data_filt[data_filt$down,]) >=1) {  # no highlighting if no down regulated genes
          points(data_filt[data_filt$down, column_name_fc], -log10(data_filt[data_filt$down, column_name_p]), 
                 pch=21, bg=color_down_reg, col="black", cex=0.7)
          mtext(text=c(expression(symbol("\257")), paste("  ", nrow(data_filt[data_filt$down,]))), side= 3, adj =0, col = color_down_reg)
        } 
 
       if(nrow(data_filt[data_filt$up,]) >=1) {  # no highlighting if no up regulated genes
         points(data_filt[data_filt$up, column_name_fc], -log10(data_filt[data_filt$up, column_name_p]), 
                pch=21, bg=color_up_reg, col="black", cex=0.7) 
         mtext(text=c(expression(symbol("\255")), paste(nrow(data_filt[data_filt$up,]), "  ")), side= 3, adj =1, col = color_up_reg)
       } 
  }

  return(data_filt)
} # end function definition


