args <- commandArgs(T)
# help file
if(length(args)<5){
  stop("Usage:
       Rscript sufficiency_test.R input threshold, num, objective, y-label
       input: the path of input file, combined randomsparse table
       threshold: percentage threshold to the objective's maximal values
       num: the number of repeats of the experiments
       objective: label showing the objective being optimized
       label: 'number of genes' or 'number of reactions'")
}


input <- args[1]
threshold <- args[2]
num <- as.numeric(args[3])
label <- args[4]
label_y <- args[5]

check_consistency <-function(now, prev, name, step){
  if(length(setdiff(prev, now))>0){
    print(paste("Previous", name, "reaction has disappeared at:", step))
  }
}

library(ggplot2)

sufficiency_test <- function(input,threshold,num,label){
  biomass <- read.table(input,sep="\t", as.is=T, row.names = 1)
  count <- matrix(nrow=3, ncol=num-1)
  core_prev <- c()
  flex_prev <- c()
  zero_prev <- c()
  for (i in seq(2, num)){
    #index <- sample(1:num, i)
    index <- 1:i
    cond <- rowSums(biomass[,index])
    core <- names(which(cond==i))  # list
    flex <- names(which(cond<i & cond>0)) # list
    zero <- names(which(cond==0))
    check_consistency(core_prev, core, "core", i)
    check_consistency(flex, flex_prev, "flex", i)
    check_consistency(zero_prev, zero, "zero", i)
    core_prev <- union(core, core_prev)
    flex_prev <- union(flex_prev, flex)
    zero_prev <- union(zero_prev, zero)
    # dir.create(paste0("sufficiency_",threshold))
    list_suff <- list(core=core, flex=flex, zero=zero)
    count[,i-1] <- sapply(list_suff,length)
  }
  rownames(count) <- c("core", "flex", "zero")
  colnames(count) <- 2:num
  print(colSums(count))
  write.table(count, paste0(input,"_",threshold,"_",num,"_sufficient_test.tsv"), sep="\t", quote=F, col.names = NA)

  suff <- data.frame(y = c(count["core",], count["flex",], count["zero",]), x = rep(2:num,3), group = c(rep("Core-essential", num-1), rep("Flexible",num-1), rep("Non-essential",num-1)))

  pdf(paste0(input,"_",threshold,"_",num,"_sufficient_test.pdf"), width=6, height=5)
  p <- ggplot(suff, aes(x=x, y=y, color=group)) +
    geom_point(size=1, alpha=0.6)
    # geom_line()
  p <- p +
    ggtitle(paste("Viability threshold of maximum", label, threshold,"%")) +
    xlab("Number of experiments") +
    ylab(label_y) +
    ylim(0,400) +
    theme_bw() +
    theme(title =element_text(size=12),
          axis.text.x = element_text(color = "grey20", size = 12),
          axis.text.y = element_text(color = "grey20", size = 12),
          axis.title.x = element_text(color = "grey20", size = 12),
          axis.title.y = element_text(color = "grey20", size = 12),
          legend.title=element_blank(),
          legend.text=element_text(size=10))
  plot(p)
  dev.off()

}

 # suff_99.99 <- sufficiency_test(threshold = 99.99, num = 1000)
 # suff_1 <- sufficiency_test(1, 1000)
# suff_99.99_3000 <- sufficiency_test(threshold = 99.99, num = 3000)
# suff_1_3000 <- sufficiency_test(threshold = 1, num = 3000)
#
# for(n in (seq(10,90,length.out=9))){
#   print(n)
#   suff_test <- sufficiency_test(threshold = n, num = 3000)
# }


# suff_99.99_1000 <- sufficiency_test(input = "../../../minimal_network/ethanol_strain/E1-mn-99.99_combine_1000_NOtfba.tsv", threshold = 99.99, num = 1000, "Ethanol production")


sufficiency_test(input, threshold, num, "Ethanol production")
