NUMTREES<-1
thetime<-as.numeric(Sys.time())
n<-16
curtree<-1
PTHRESH_STR<-c("0","1e-4","5e-4",".001",".005",".01",".05")
PTHRESH<-c(0,.0001,.0005,.001,.005,.01,.05)
options(shiny.maxRequestSize=100*1024^2)
genebardata <- read.csv("genebardata.csv")
phenonames<-as.character(genebardata$phenovar)

oldpsx=matrix(rep(list()),16,nrow=16,ncol=1)
oldpsy=matrix(rep(list()),64,nrow=16,ncol=4)

selrow1=list()
