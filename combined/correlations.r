####### Correlations between effect sizes of CpGw with lowest p-value (top 100) in different summary statistics
# M.Suderman

# load summary statistics (CpG, estimate, p-value)
# create list of analysis with lists of CpGs available in all analysis (ewas.list)

get.top.sites <- function(ewas,n) 
  rownames(ewas)[order(ewas$p.value)[1:n]]

top.sites <- sapply(ewas.list, get.top.sites, n=100)

top.corr <- sapply(1:length(ewas.list), 
                   function(i) 
                     sapply(1:length(ewas.list), 
                            function(j) 
                              cor(ewas.list[[i]][top.sites[,i], "estimate"],
                                  ewas.list[[j]][top.sites[,i], "estimate"], 
                                  use="p")))
print(top.corr)
