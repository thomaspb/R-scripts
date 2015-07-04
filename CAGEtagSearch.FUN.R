# CAGEtagSearch.FUN.R



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Progress bar settings
total <- nrow(CAGEinputX) 
pb <- txtProgressBar(min = 0, max = total, style = 3)

remove(CAGE.result)
CAGE.result = list()
for (i in 1:nrow(CAGEinputX)) 
{
  ## Decompose the CAGE tag into its components
  x <- as.character(rownames(CAGEinputX[i,]))
  m <- regexec("^(([^:]+):)?([0-9]+))?(..([0-9]+))?(,-|+.*)", x)
#  m
#  regmatches(x, m)
  #        
  CAGEtag_parts <- function(x) {
    m <- regexec("^(([^:]+):)?([0-9]+))?(..([0-9]+))?(,-|+.*)", x)
    parts <- do.call(rbind, lapply(regmatches(x, m), `[`, c(3L, 4L, 6L, 7L)))
    colnames(parts) <- c("chromosome","start","stop","strand")
    parts}
  #        
  CAGEtaginfo <- CAGEtag_parts(x)
  # head(CAGEtaginfo)
  # if (CAGEtaginfo[,1] == "chr17" & CAGEtaginfo[,2] >= 10009524 & CAGEtaginfo[,3] <= 10009534 & CAGEtaginfo[,4] == ",-" ) 
  if (CAGEtaginfo[,1] == TargetChromosome & CAGEtaginfo[,2] >= TargetStart & CAGEtaginfo[,3] <= TargetStop & CAGEtaginfo[,4] == TargetStrand ) 
  { 
    # print(CAGEinputX[i,]) 
    CAGE.result[[i]] = CAGEinputX[i,]    
    # update progress bar
    
  }  # END OF IF STATEMENT
  #
  setTxtProgressBar(pb, i)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 