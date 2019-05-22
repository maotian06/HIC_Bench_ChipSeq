#$ -S /usr/bin/Rscript

# libraries
library(optparse)
library(reshape)
library(preprocessCore)

## ##############################################
## run
## ##############################################

run = function(max_dist,pseudo_bind,pseudo_expr,qnorm_bind,qnorm_expr,file='') {

  # compute fold-changes on binding data
  i1 = 12
  i2 = which(colnames(dbind)=='feature')-1
  fcbind_func = function(i) { log2((dbind[,i+1]+pseudo_bind)/(dbind[,i]+pseudo_bind)) }               # FC function, with pseudocount
#  if (qnorm_bind==TRUE) dbind[,i1:i2] = normalize.quantiles(as.matrix(dbind[,i1:i2]))                # qnorm BEFORE fold-change computation (NOT USED)
  fc_bind = fcbind_func(seq(i1,i2,by=2))																															# FC computation: R vs D
  if (qnorm_bind==TRUE) fc_bind[,1:ncol(fc_bind)] = normalize.quantiles(as.matrix(fc_bind))						# qnorm on fold-changes
  colnames(fc_bind) = gsub('[.]R[.].*$','',colnames(fc_bind))		       																# fix sample names
  p = ((dbind$insideFeature=='overlapStart')|((dbind$insideFeature=='upstream')&(abs(dbind$shortestDistance)<=max_dist))|(dbind$insideFeature=='includeFeature'))     # select features
  fc_bind = melt(cbind(dbind[p,i2+2],fc_bind[p,]))                                                    # melt
  colnames(fc_bind)[1:2] = c('GENE','PATIENT')                                                        # fix column labels

  # compute fold-changes on expression data
  fcexpr_func = function(i) { log2((expr[,i+1]+pseudo_expr)/(expr[,i]+pseudo_expr)) }                 # FC function, with pseudocount
  if (qnorm_expr==TRUE) expr[,2:ncol(expr)] = normalize.quantiles(as.matrix(expr[,2:ncol(expr)]))     # qnorm BEFORE fold-change computation
  fc_expr = fcexpr_func(seq(2,ncol(expr)-1,by=2))                                                     # FC computation: R vs D
#  if (qnorm_expr==TRUE) fc_expr[,1:ncol(fc_expr)] = normalize.quantiles(as.matrix(fc_expr))          # qnorm AFTER fold-change computation (NOT USED)
  colnames(fc_expr) = gsub('[._]R$','',colnames(fc_expr))                                             # fix sample names
  fc_expr = melt(cbind(expr[,1],fc_expr))                                                             # melt
  colnames(fc_expr)[1:2] = c('GENE','PATIENT')                                                        # fix column labels

  # merge the two datasets
  merged_data = merge(fc_bind,fc_expr,by=c('GENE','PATIENT'))                                         # merge tables, by GENE/PATIENT
  colnames(merged_data)[3:4] = c('BIND','EXPR')                                                       # fix column labels

  # save results
  if (file!="") save(merged_data,file=file)
  
  # evaluate
  z = (abs(merged_data$BIND)>2.0)&(abs(merged_data$EXPR)>2.0)                                         # keep extreme folc-changes only
  #sum(z)
  return(cor(merged_data$BIND[z],merged_data$EXPR[z]))                                                # return correlation (across all samples)
}


run_parameter_exploration = function()
{
  for (max_dist in c(1000,5000))
    for (pseudo_bind in c(0.1,1.0))
      for (pseudo_expr in c(10.0,50.0))
        for (qnorm_bind in c(TRUE,FALSE))
          for (qnorm_expr in c(TRUE,FALSE)) {
            corr = run(max_dist,pseudo_bind,pseudo_expr,qnorm_bind,qnorm_expr)
            cat(paste(max_dist,pseudo_bind,pseudo_expr,qnorm_bind,qnorm_expr,corr,sep='\t'))
            cat('\n')
          }
}


## ##############################################
## MAIN
## ##############################################

# command-line
args <- commandArgs(trailingOnly=T)

# process command-line arguments
option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-o","--output-file"), default="", help="Output RData file (required) [default=\"%default\"].")
)
usage = 'diffbind-integrate [OPTIONS] DIFFBIND-MATRIX COUNT-MATRIX';

# get command line options (if help option encountered print help and exit)
arguments = parse_args(args=args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf))
opt = arguments$options
fout = opt$'output-file'
files = arguments$args
if (length(files)!=2) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
file_diffbind = files[1]
file_counts = files[2]

# read matrices
if (opt$verbose) write("Loading tables...",stderr())
dbind = read.csv(file_diffbind)
expr = read.table(file_counts,header=T,sep='\t',check.names=F)
#expr = read.csv(file_counts)

# merge diffbind and expression tables
if (opt$verbose) write("Merging tables...",stderr())
X = merge(y=dbind,x=expr,by.y='external_gene_name',by.x='Gene',all.y=TRUE)
write.csv(X,file=fout,row.names=F)

# done
if (opt$verbose) write("Done.",stderr())
quit(save='no')



# parameter exploration
run_parameter_exploration()

# output file
fout.rdata = paste(opt$'output-file','.RData',sep='')
fout.pdf = paste(opt$'output-file','.pdf',sep='')

# use "optimal" parameters
print(run(max_dist=1000,pseudo_bind=0.1,pseudo_expr=10,qnorm_bind=TRUE,qnorm_expr=TRUE,file=fout.rdata))
#print(run(max_dist=1000,pseudo_bind=0.1,pseudo_expr=10,qnorm_bind=TRUE,qnorm_expr=FALSE,file=fout.rdata))

# load, process and plot
load(fout.rdata)
fplot = function(merged_data,patient,cutoff_bind,cutoff_expr,N)
{
  z = (merged_data$PATIENT==patient)&(abs(merged_data$BIND)>=cutoff_bind)
  m = merged_data[z,]
  N = min(N,nrow(m))
  expr_pos_cutoff = max(+cutoff_expr,m$EXPR[order(m$EXPR,decreasing=TRUE)[N]])
  expr_neg_cutoff = min(-cutoff_expr,m$EXPR[order(m$EXPR,decreasing=FALSE)[N]])
  boxplot(m[m$EXPR>=expr_pos_cutoff,]$BIND,m[m$EXPR<=expr_neg_cutoff,]$BIND,main=patient,names=c('UP','DN'),ylab='DiffExpr',cex=0.5)
#  boxplot(m[m$BIND>=cutoff_bind,]$EXPR,m[m$BIND<=-cutoff_bind,]$EXPR,main=patient,names=c('UP','DN'),ylab='DiffBind')
}

pdf(fout.pdf)
for (cutoff_expr in c(0.01,1.0,1.5,2.0)) {
  for (cutoff_bind in c(1.0,1.5,2.0)) {
    for (N in c(100,200,500,1000,Inf)) {
      par(mfrow=c(2,5),mar=c(2,2,4,2),oma=c(1.5,2,1,1))
      for (p in as.vector(unique(merged_data$PATIENT))) fplot(merged_data,p,cutoff_bind,cutoff_expr,N=N)
      mtext(paste("e-cutoff=",cutoff_expr,", b-cutoff=",cutoff_bind,", N=",N,sep=''),outer=TRUE,cex=1,line=-0.5)
    }
  }
}
dev.off()




