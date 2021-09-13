args <- commandArgs(T)
# help file
if(length(args)<5){
  stop("Usage:
       Rscript robust-fva-parsing_modified.R infile varylabel
       infile: the path of input file
       varylabel: the x axis lable, based on the optimized reaction, e.g. if vary etahnol production, 
              the varylabel could be 'Ethanol production'.
       PlotRxns: a tab-separated two-column input file that specify the reactions need to be ploted, 
              1st column is reaction ID (e.g. R01196), 2nd column is reaction label (e.g.POR)
       PanelDimen: the number of panels, c(row, col), in the out .png file. For example, 
              '3,2' means there will be 3 plots per row and two columns in final .png.
       outname: suffix of output image file. If 'infile' name is 'E1-robustness-fva-result.tsv' and 
              '-plot' is given as outname, then the final image output name is 
              'E1-robustness-fva-result.tsv-plot.png'. By default, outname is '-plot-rxn'.
       "
       
      )
}


infile <- args[1]
varylabel <- args[2]
PlotRxns <- args[3]
PanelDimen <- args[4]

# set suffix of output image file name
if (is.na(args[5])) {
  outname = '-plot-rxn'
} else {outname = args[5]}
outname

#rm(list=ls())
########## Function ########################################
#-function for producing plots based on the robustness-fva analysis
makeplot<-function(outfile,paneldim,steps,varylabel,rxnmap,r) {
  #-outfile: the output .png file
  #-paneldim: the number panels, c(row, col) in the out .png file
  #-steps: the number of steps in a simulation
  #-varylabel: the label of the varying axis
  #-rxns: the list of reactions to plot
  #-r: read.table import from the original output of robust-fva.tsv.py
  png(outfile,res=300,width=700*paneldim[2],height=700*paneldim[1])
  # pdf(outfile,res=300,width=700*paneldim[2],height=700*paneldim[1])
  par(mfrow=paneldim)
  for(i in 1:dim(rxnmap)[1]){
    rxn = as.matrix(rxnmap[i,1])
    if(! rxn %in% r[,1]){
      plot.new()
      next
    }
    header = as.matrix(rxnmap[i,2])
    fluxes = matrix(0,nrow=steps,ncol=dim(r)[2]-1)
    ssize = length(unique(r[,1]))
    for(j in 0:(steps-1)){
      start = j*ssize+1
      end = (j+1)*ssize
      tmp = as.matrix(r[start:end,2:dim(r)[2]])
      rownames(tmp) = r[start:end,1]
      fluxes[j+1,] = tmp[rxn,]
    }
    x = fluxes[,1]
    y1 = fluxes[,2] #min
    y2 = fluxes[,3] #max
    # plotmin = floor(min(y1,y2)/10) * 10
    # plotmax = ceiling(max(y1,y2)/10) * 10
    plotmin = (floor(min(y1,y2)*100)/100)
    plotmax = (ceiling(max(y1,y2)*100)/100)
    plotbound_pre = max(abs(plotmin), abs(plotmax))
    plotbound = max(abs(plotmin), abs(plotmax)) * (1+0.3)
    x_bound = max(x) * (1+0.22)               
    plot(x, y1, col='blue', cex=0.5, main=header, xlab=varylabel, ylab='Flux Range',xlim=c(0, x_bound), ylim=c(-plotbound,plotbound))
    rows = nrow(fluxes) # test how many rows the fluxes has
    start_upper = round(fluxes[1, 3], 2)
    start_lower = round(fluxes[1, 2], 2)
    end_upper = round(fluxes[rows, 3], 2)
    end_lower = round(fluxes[rows, 2], 2)
    # upper_bound = paste("(", round(fluxes[1, 3], 2), ", ", round(fluxes[rows, 3], 2), ")", sep="") # upper bound at start and end point
    # lower_bound = paste("(", round(fluxes[1, 2], 2), ", ", round(fluxes[rows, 2], 2), ")", sep="") # lower bound at start and edn point
    # text(x=0, y=max(y2), paste("upper bound (start, end):", upper_bound, sep = " "), cex=0.65, adj=c(0,0))
    # text(x=0, y=min(y1), paste("lower bound (start, end):", lower_bound, sep = " "), cex=0.65, adj=c(0,1))
    text(x=x[1], y=y2[1] + plotbound_pre * 0.12, labels=start_upper, cex=0.8, font=2, adj=c(0,0))
    text(x=x[1], y=y1[1] - plotbound_pre * 0.12, labels=start_lower, cex=0.8, font=2, adj=c(0,1), col="blue")
    text(x=tail(x, n=1), y=tail(y2, n=1) + plotbound_pre * 0.15, labels=end_upper, cex=0.8, font=2, adj=c(0,0))
    text(x=tail(x, n=1), y=tail(y1, n=1) - plotbound_pre * 0.20, labels=end_lower, cex=0.8, font=2, adj=c(0,1), col="blue")
    polygon(c(x,rev(x)), c(y2,rev(y1)), col='grey')
    # lines(x, y2, col='red', cex=1, lty=1)
    abline(h=0, col = 'black', lty=2)
  }
  dev.off()
}


########## Main Script ########################################
#infile = '../../../robustness-fva/result/3-dGOR-robustness-fva_vary_sum-products_opt_biomass.tsv'
r = read.table(infile,sep="\t")
steps = length(r[,1])/length(unique(r[,1])) #the number of steps in a robustness analysis
#varylabel = 'Products Production (mM)'

# #-carbon plot
# outfileC = paste(infile,'-Carbon.png',sep='')
# carbonrxns = c('sink_biomass', 'TP_CO2', 'Lac-symp', 'PYRt2', 'ACt6')
# carbonlabels = c('Biomass', 'CO2', 'Lactate', 'Pyruvate', 'Acetate')
# carbonrxnmap = data.frame(carbonrxns,carbonlabels)
# makeplot(outfileC,c(2,3),steps,varylabel,carbonrxnmap,r)
# 
# 
# #-redox plot
# outfileH2 = paste(infile,'-H2.png',sep='')
# h2rxns = c('BF-Nfn', 'BF-H2ase', 'MBH', 'R01196','R01061','R07159', 'TP_H2', 'EX_C00080[e]', 'Rnf')
# H2labels = c('BF-Nfn', 'BF-H2ase', 'MBH', 'POR','GAPDH','GOR', 'H2', 'H+', 'Rnf')
# h2rxnmap = data.frame(h2rxns,H2labels)
# makeplot(outfileH2,c(3,3),steps,varylabel,h2rxnmap,r)

# carbon+redox plot based on rxns listed in argument "PlotRxns"
outfileC = paste(infile, outname, ".png", sep='')
# rxns = read.table(PlotRxns, sep="\t") # 1st column is rxn ID, 2nd column is rxn label
# rxns = strsplit(target, ",")[[1]]
# labels = strsplit(label, ",")[[1]]
rxnmap = read.table(PlotRxns, sep="\t")
makeplot(outfileC,as.numeric(strsplit(PanelDimen, ",")[[1]]),steps,varylabel,rxnmap,r)




