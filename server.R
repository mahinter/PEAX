library(shiny)
library(DT)
library(parallel)
library(miscTools)
library(ggplot2)
#library(multiMiR)
source('df2html.R')

curboxplot=0
oldselmir=0
oldselmir1=0 
oldselmir2=0 
oldsearch=0 
oldsearchall=0 
olddoofat=rep(0,256) 
oldofat=data.frame()
lastofat=0 
oldfindbest=0 
oldloaddemo=0 
loadtree=0
demoloads=0 
numtests=0
curfilt=0
selopts=list()
ppvals=data.frame()
ppsize=data.frame()
normals=list()
defs=NULL

mcadply <- function(X, FUN, ...) {
    # Runs multicore lapply with progress indicator and transformation to
    # data.table output. Arguments mirror those passed to lapply.
    #
    # Args:
    # X:   Vector.
    # FUN: Function to apply to each value of X. Note this is transformed to 
    #      a data.frame return if necessary.
    # ...: Other arguments passed to mclapply.
    #
    # Returns:
    #   data.table stack of each mclapply return value
    #
    # Progress bar code based on http://stackoverflow.com/a/10993589
    require(parallel)
    require(plyr)
    require(data.table)
    
    local({
        f <- fifo(tempfile(), open="w+b", blocking=T)
        if (inherits(parallel:::mcfork(), "masterProcess")) {
            # Child
            progress <- 0
            print.progress <- 0
            while (progress < 1 && !isIncomplete(f)) {
                msg <- readBin(f, "double")
                progress <- progress + as.numeric(msg)
                # Print every 1%
                if(progress >= print.progress + 0.01) {
                    cat(sprintf("Progress: %.0f%%\n", progress * 100))
                    print.progress <- floor(progress * 100) / 100
                }
            }
            exit()
        }
        
        newFun <- function(...) {
            writeBin(1 / length(X), f)
            return(as.data.frame(FUN(...)))
        }
        
        result <- as.data.table(rbind.fill(mclapply(X, newFun, ...)))
        close(f)
        cat("Done\n")
        return(result)
    })
}

rep.row<-function(x,n){
     matrix(rep(x,each=n),nrow=n)
}

get_colors<-function() {
  return (c("red","blue","green","gray","orange","yellow","purple"))
}
shinyServer(function(input, output) {
   output$tablespacing <- renderText({
    "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n."
   })
   output$tablespacing2 <- renderText({
    "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n."
   })

   testcor2<-reactive({
     return(testcor(2))
   })

   testcor1<-reactive({
     return(testcor(1))
   })

    datasetInput<-reactive({
     inFile <- input$file1

     if (!is.null(inFile)) {
     #message("READING DATASET") 
     x<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
     }
     else if (input$loaddemo > 0) {
        #message("DEFAULT DATASET") 
        x<-read.csv("./data/peax_demo_pheno.csv", header=input$header, sep=input$sep, quote=input$quote)
     }
     else
     {  
	return(NULL) }
     return(x)}) 

   mirInput<-reactive({
     #input$search
     #message("READING MIRINPUT") 
     mirFile <- input$mirfile

     if (!is.null(mirFile))
       y<-read.csv(mirFile$datapath,check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
     else if (input$loaddemo > 0) {
        #message("DEFAULT DATASET") 
        y<-read.csv("./data/peax_demo_mrna.csv", header=TRUE, sep=",", quote="\"")
     }
     else 
       return(NULL)
     
     if (nrow(y)>50000) {
 	y<-y[1:50000,]
     }
     y1<-y
     if (input$transpose) {
       message("Transposing Exp1...")
       yn=y[,1]
       cn=colnames(y1)[-1]
       y<-as.data.frame(t(y[,2:ncol(y)]))
	y<-setNames(y,yn)
	y<-cbind(cn,y)
	colnames(y)[1]<-"Pat"
     }
     #if (input$search > oldsearch)
     #{
#	#message("New Search!")
  #   }
   #  oldsearch<<-input$search
     return(y)}) 
   evidInput<-reactive({
     #message("READING EVIDENCE") 
     evidFile <- input$evidence

     if (is.null(evidFile))
       return(NULL)
     z<-read.csv(mrnaFile$datapath, check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
     return(z)})


   defInput<-reactive({
     defFile <- input$deffile

     if (!is.null(defFile))
	z<-read.csv(defFile$datapath, check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
     else if (input$loaddemo > 0) {
 	#message("DEFAULT DATASET") 
        z<-read.csv("./data/peax_demo_def.csv", header=TRUE, sep=",", quote="\"")
     }
     else
       return(NULL)
     return(z)
   })

   mrnaInput<-reactive({
     #message("READING MRNAINPUT") 
     mrnaFile <- input$mrnafile

     if (!is.null(mrnaFile))
	z<-read.csv(mrnaFile$datapath, check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
     else if (input$loaddemo > 0) {
 	#message("DEFAULT DATASET") 
        z<-read.csv("./data/peax_demo_mir.csv", header=TRUE, sep=",", quote="\"")
     }
     else
       return(NULL)
     if (input$transpose) {
       message("Transposing Exp2...")
       zn=z[,1]
       z<-as.data.frame(t(z[,2:ncol(z)]))
 	z<-setNames(z,zn)
     }
     return(z)}) 

rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}


################

get_selection<-function()
{
  if (length(input$testTbl1_1_rows_selected)>0)
	{
	 selrow1<<-as.numeric(input$testTbl1_1_rows_selected)
	 return (selrow1)
        }
  if (length(selrow1) > 0)
	return (selrow1)
  ppvals<<-data.frame()
  ppsize<<-data.frame()
  return (list())
}

get_selection_names<-function()
{
    sels<-get_selection()
    selnames<-colnames(mirInput())[as.numeric(sels)+1] #add 1 to skip first column
    retnames<-data.frame(sels=sels,names=selnames)
    return(retnames)
}

get_splitpheno<-function(x,node,tr,cols,vals)
{
  message(paste0("DEBUG:get_splitpheno",node," tree:",tr)) 
  retpheno<-x
  #
#  1
#2  3
#45  67
  if ((cols==0) || (length(cols)==0)) { 
    curvar<-paste0("var",tr,"_",node%/%2)
    cursld<-paste0("range_slider",tr,"_",node%/%2)
    if (is.null(input[[curvar]]) || is.null(input[[cursld]])) {
 	#message("NULL HERE")
	return(x)
     } 
    slval<-input[[cursld]]
    curvar<-input[[curvar]] 
  }
  else {
    curvar<-cols[node%/%2]
    slval<-vals[node%/%2]
  }
  message(paste0("curvar:",curvar))
  if (node==1) {
	return(x)
  } 
  #x
 
  if (node>1) {
    retpheno<-get_splitpheno(x,node%/%2,tr,cols,vals) #get parent pheno
    if (is.null(retpheno) || nrow(retpheno)==0) {
      retpheno<-retpheno
    }
    else if (node%%2==0) { #odd/left branch
      #message("Debug left")
      retpheno<-retpheno[retpheno[curvar]<slval,]
      #message("Debug leftdone")
    }
    else {
      #message("Debug right")
      retpheno<-retpheno[retpheno[curvar]>=slval,]
      #message("Debug rtdone")
    }
  }
  else
  {
   retpheno<-x  ## everybody
  }
  ds<-input[[paste0("depth_slider",tr)]]
#  rr$Group=4
  # if it's a leaf, add the group annotation
  leaf_thresh=2^ds
  
  if (!is.null(retpheno) && nrow(retpheno)>0) {
    retpheno$Group=node%%leaf_thresh+1
  }
  else {
   message("DEBUG:returning NULL retpheno")
  }
  return(retpheno)
}

getPhenoExp<-function(x, curtr, curboxplot) {

y <- mirInput()

ds1<-input$depth_slider1
ds2<-input$depth_slider2

if (curtr==1) {
  oldselmir=oldselmir1
  ds<-input$depth_slider1
  tbl<-input$testTbl1_1
}
else
{
  oldselmir=oldselmir2
  ds<-input$depth_slider2
  tbl<-input$testTbl2_1
}

 #selmir<-tbl[curtr]
 #curboxplot<-input$testTbl1_1_row_last_clicked
 
 selmir<-curboxplot
 # BROKEN here? selmir<-"rand.miR.101"
 #message("getPhenoExp")
 #message(selmir)
 #message(y)
 #message("getPhenoExp2")
 if (is.null(selmir) || is.na(selmir)) {
   if (is.null(oldselmir) || oldselmir==0)
     selmir = colnames(y)[2]
   else
     selmir = oldselmir
 }
 else {
 	selmir<-colnames(y)[as.numeric(selmir)+1] # convert string to numeric, and add 1 for first col
 	#selmir=selmir
 } 

 #make a simple pheno vector with PT column and pheno value
 themir<-cbind(y[1],y[selmir])
 colnames(themir)[2]<-"Pheno"

 #message("DEBUG:getPheno3d")
 x<-merge(x,themir,all.x=TRUE,by=1) #merge pheno by Pt # in first col
 message("DEBUG:getPhenoExp3d")
 return(x)
}



getPheno<-function(curtr) {
message(paste0("DEBUG:getPheno tree:",curtr))
for (t in 1:2) {
for (i in 1:8) {
  a<-input[[paste0("range_slider",t,"_",i)]]
  #message(a)
}}
ds1<-input$depth_slider1
ds2<-input$depth_slider2

if (curtr==1) {
  oldselmir=oldselmir1
  ds<-input$depth_slider1
  tbl<-input$testTbl1_1
}
else
{
  oldselmir=oldselmir2
  ds<-input$depth_slider2
  tbl<-input$testTbl2_1
}

if (is.null(ds)) {
  return(NULL)
}

x <- datasetInput()
#message("DEBUG:getPheno2")
if (is.null(x)) {
  #message("DEBUG:datasetInput is null")
  return(NULL)
}
x<-x[c(colnames(x)[1],sapply(1:(2^ds-1),function(x) input[[paste0("var",curtr,"_",x)]]))]
x <- x[complete.cases(x), ]
#message("DEBUG:getPheno3")
y <- mirInput()

# bind all leaves
#e.g. depth 3
# 8 leaf nodes, labeled 8 through 15

leaf_nodes<-2^(ds)
ret<-NULL
#message("DEBUG:getPheno3e")
for (lf in leaf_nodes:(leaf_nodes*2-1)) {
  ret<-rbind(ret,get_splitpheno(x,lf,curtr,0,0))
}
#message("DEBUG:getPheno3f")
if (!is.null(ret)) {
  ret<-ret[with(ret,order(ret[,1])),]
}
#message("DEBUG:getPheno3g")
#message("DEBUG:Leaving getPheno")
return(ret)
}


get_pvals<-function(x,y,ds,fact,cols,vals,genes)
{
leaf_nodes<-2^(ds)
ret<-NULL
#message("DEBUG:getPheno3e")
message(paste0("cols:",cols))
message(paste0("vals:",vals))
for (lf in leaf_nodes:(leaf_nodes*2-1)) {
  ret<-rbind(ret,get_splitpheno(x,lf,tr,cols,vals))
}

if (!is.null(ret)) {
  ret<-ret[with(ret,order(ret[,1])),]
}

print(ret)
 message("DEBUG:HERE2")
y<-cbind(y,ret["Group"])
y$Group<-factor(y$Group)
ps<-data.frame()
groupsize<-min(data.frame(table(y$Group))$Freq)
for (gene in genes) {
  form<-as.formula(paste(gene,"Group",sep="~"))
  print(form)
  if (nlevels(y$Group)<2)
	curp<-1
  else {
     if (input$testtype=='Kruskal-Wallis')
       curp<-kruskal.test(form,y)$p.value
     else
       curp<-anova(lm(form,y))$"Pr(>F)"[1]
  }
  message(paste0("curp:",curp))
  message(paste0("gene:",gene))
  message(paste0("vals:",vals))
  message(y)

  ps<-c(ps,curp)
}

ps<-data.frame(t(ps))
colnames(ps)<-genes
rownames(ps)<-as.numeric(vals[fact])
ps<-cbind(ps,rep(groupsize,nrow(ps)))
colnames(ps)[ncol(ps)]<-"minsize"
return (ps)
}

get_ofat<-function(x,y,ds,fact,val,genes,doall)
{
#temporarily fix this to median
#val<-colMedians(x[fact])


# find low and hi vals based on median.  can generate more
message(paste0("fact:",fact))
ecf<-ecdf(x[,fact])
pct<-ecf(val)
low<-pct*0.8
hi<-pct+(1-pct)*0.2

low<-quantile(x[,fact],low)
hi<-quantile(x[,fact],hi)


cols<-colnames(x)[2:length(colnames(x))] #all geno columns, so skip first one

# temporarily set to median, instead of input vals
vals<-(apply(x[,cols,drop=F],2,median))
vals[]<-sapply(1:(2^ds-1),function(x) input[[paste0("range_slider",1,"_",x)]])
vals[fact]<-val
#allvals<-sort(x[,fact])
allvals<-x[,fact]
message(paste0("Starting all vals!",allvals))
message(paste0("cols!",cols))
message(paste0("colnames!",colnames(x)))
message(paste0("vals!",vals))

of<-data.frame()

if (doall) {
  for (val in unique(unlist(allvals, use.names=FALSE))) {
    message("************************Doing value:",val)
    vals[fact]<-val
    curof<-get_pvals(x,y,ds,fact,cols,vals,genes)
    of<-rbind(of,curof)
  }
}
else {
#for (val in allvals) {
  message("Doing value:",val)
of<-get_pvals(x,y,ds,fact,cols,vals,genes)
#vals[fact]<-low
#curof<-get_pvals(x,y,ds,fact,cols,vals,genes)
#of<-rbind(of,curof)
#vals[fact]<-hi
#curof<-get_pvals(x,y,ds,fact,cols,vals,genes)
#of<-rbind(of,curof)
}

#}
message("ending all vals!")

return (of)
}



  getPheno1<-reactive({
    return(getPheno(1))
  })
  
  getPheno2<-reactive({
    return(getPheno(2))
  })



   testcor <-function(curtr)
   {
   #if (tr==1) {
   #  data<-getPheno1()
   #}
   #else {
   #  data<-getPheno2()
   #}
   #if (is.null(data)) {
   #	return(NULL)
   #}
   #message("DEBUG:Enter testcor()...")
   z<-mrnaInput()
   y<-mirInput()
   
   if (is.null(z) || is.null(y)) {
 	return(NULL)
   }
   y<-y[y[,1] %in% z[,1],]
   z<-z[z[,1] %in% y[,1],]
   
   # do we need Group for corr?
   #z$Group<-data$Group
 
   # remove patients (col 1) before calculating
   
 if (curtr==1)
   tbl<-input$testTbl1_1
 else
   tbl<-input$testTbl2_1
 
  selmir<-tbl[curtr]

   if (is.null(selmir) || is.na(selmir)) 
   { 
	#message("selmir undefined in testcor")
   	selmir<-oldselmir
   }
   #message(paste0("DEBUG:Starting correlation...:",selmir))
   pvals<-cor(y[selmir],z[-1],method="spearman")
   #message("DEBUG:Ending correlation...")
   #pvals<-pvals[,abs(pvals)>0.6]
   pvals<-sort(pvals[1,])
   #nm<-names(pvals)
   nm<-names(pvals)
   zcor=data.frame(Probe = nm,Spearman=round(as.numeric(pvals),4))
   #message("DEBUG:Leaving testcor()...")
   return(zcor)
   }         
  

   testdf1 <- reactive(
   {
     return(testdf(1))
   })
   testdf2 <- reactive(
   {
     return(testdf(2))
   })

   testdf<-function(tr) {
   ts<-input$tabSelected
   #message(ts)
   a <- input$range_slider1_1

   if (tr==1) {
     data<-getPheno1()
   }
   else {
     data<-getPheno2()
   }
   if (is.null(data)) {
        #message("DEBUG:NULL data in getPheno")
 	return(NULL)
   }
   #message(paste0("DEBUG:Enter testdf()... with tree:",tr))
   y<-mirInput()
   data<-data[data[,1] %in% y[,1],]
   y<-y[y[,1] %in% data[,1],]
 
   y$Group<-data$Group
   y$Group<-as.factor(y$Group)
   yaov<-0
   
   #message("DEBUG:begin Anova step")

   y$Group<-data$Group
   y<-data.frame(as.matrix(y[-1]))
   y$Group<-as.factor(y$Group)
   Klist<-sort(unique(y$Group))
   #message(paste0("DEBUG:Klist=",Klist))
   K<-length(Klist)
   ### FILTER HERE
   #Yib<-mclapply(Klist,function(tmp) colMedians(y[y$Group==tmp,][-ncol(y)]))
   
   if (input$searchall > oldsearchall) {
	oldsearchall<<-input$searchall
	oldsearch<<-input$search
	curfilt<<-NULL  
	selrow1<<-list()
 	oldofat<<-data.frame()
   }
 
   if (input$findbest > oldfindbest) {
	oldfindbest<<-input$findbest
   }
   
   if (input$search > oldsearch)
   {
	#message(paste0("Search:",input$search," ",oldsearch))
	oldsearch<<-input$search
        curfilt<<-input$testTbl1_1
   	if (is.null(curfilt)) {
     	  curfilt<<-oldselmir1
   	}
	#message(paste0("New Search!",curfilt))
   }
   	if (!is.null(curfilt)) {
   		#Ydf<-data.frame(Yib) 
		#curthresh<-input$fold_thresh/100
   		#Yfilt<-mclapply(Klist,function(tmp) median(y[y$Group==tmp,curfilt]))
   		#Yfilt<-(Yfilt<(1-curthresh))*-1+(Yfilt>(1+curthresh)*1)   
   		#Ydf<-data.frame((Ydf<(1-curthresh))*-1+(Ydf>(1+curthresh)*1))   
   		#found<-t(data.frame(colSums(t(Ydf)==Yfilt)))    
   		#names(found)<-colnames(y)[-ncol(y)]
   		#found<-found[found==K]
   		#found<-names(found)
   		#y<-y[found] 
   		#y$Group<-as.factor(data$Group)

	}
	else {
	}

   nm<-names(y)[-ncol(y)]
   #nm<-nm[-1]
   N<-nrow(y)
   Ybar<-colMeans(y[-ncol(y)])
   #Yib<-mclapply(Klist,function(tmp) colMeans(y[y$Group==tmp,][-ncol(y)]))
   #ni<-mclapply(Klist,function(tmp) nrow(y[y$Group==tmp,]))
   #expV<-mclapply(1:K,function(tmp) (ni[[tmp]]*(Yib[[tmp]]-Ybar)^2)/(K-1))
   expVar<-0
   #for (i in 1:length(expV)) {
   #  expVar=expVar+expV[[i]]
   #}  

   #yy<-mclapply(Klist,function(tmp) apply(y[y$Group==tmp,][-ncol(y)],2,'-',Yib[[tmp]]))
   unexpV<-0
   #for (i in 1:K) {
   #  yy<-y[y$Group==Klist[i],][-ncol(y)]-rep.row(Yib[[i]],ni[[i]])
   #  yy<-yy^2/(N-K)
   #  unexpV=unexpV+colSums(yy)
   #}  

   #Ft<-expVar/unexpV

   #anovapval<-(1-pf(Ft,df1=K-1,df2=N-K))
   message("DEBUG:end Anova step")
   
   numtests<<-numtests+length(y)
   forms <- paste(nm, " ~ Group", sep="")
   doKW<-function(eachFormula,myData) {kruskal.test(as.formula(eachFormula),data=myData)$p.value}
   doAnova<-function(eachFormula,myData) {anova(lm(as.formula(eachFormula),data=myData))$"Pr(>F)"[1]}

  if (input$testtype=="Kruskal-Wallis")
     pvals<-unlist(mclapply(forms,doKW,myData=y))
   else {
     pvals<-unlist(mclapply(forms,doAnova,myData=y))
    }
   yaov=data.frame(Probe = nm,p_val=round(as.numeric(pvals),5))
   #WHY REMOVE THE LAST ROW?yaov=yaov[-1,]
   yaov=yaov[with(yaov,order(p_val)),]
   #message("DEBUG:Leaving testdf()...")
   return(yaov)
   }         
  
    numTests <- reactive({
	a <- input$range_slider1_1
	b <- input$range_slider1_2
	c <- input$range_slider1_3
	d <- input$range_slider1_1
	e <- input$range_slider1_2
	f <- input$range_slider1_3
	g <- input$searchall
	h <- input$search
       paste0("#Tests:",numtests)
    })
    
    normalityTest <- reactive({
	a <- input$range_slider1_1
	b <- input$range_slider1_2
	c <- input$range_slider1_3
	d <- input$range_slider1_1
	e <- input$range_slider1_2
	f <- input$range_slider1_3
	g <- input$searchall
	h <- input$search
       paste0("Normality:",round(min(unlist(normals)),2))
	})

    output$numtests<-renderText({
	  paste0("  ",numTests())
	})
    output$normality<-renderText({
	  paste0("  ",normalityTest())
	})

    for (tr in 1:NUMTREES) {
    local ({
     curtr<-tr
     tbl1<-paste0("testTbl",curtr,"_1") 
     tbl2<-paste0("testTbl",curtr,"_2") 
    # messages MIR table with selectable rows and conditional formatting
    output[[tbl1]] <- DT::renderDataTable({
    tdf<-testdf(curtr)



    if (!is.null(tdf)) {
      # restrict to top 32 mirs
      tdf<-na.omit(tdf[1:64,])
      tdf[,3]<-tdf[,2]
      colnames(tdf)[3]<-colnames(tdf)[2]
      colnames(tdf)[2]<-"Link"
      colnames(tdf)[2]<-"p-val"
      tdf<-tdf[,1:2]
      mbnames<-tdf
      mbnames[,1]<-gsub("_st","",mbnames[,1])
      mbnames[,1]<-gsub("[.]","_",mbnames[,1])
      mbnames[,1]<-gsub("_star","_5p",mbnames[,1])
      #tdf[,2]<-paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",mbnames[,1]," target=_blank>",mbnames[,1],"</a>",sep="")
     
       #HTML(df2html(tdf, class = "tbl selRow", id = tbl1, 
      #             cellClass = cbind(rep(NA, nrow(tdf)), rep(NA, nrow(tdf)), ifelse(tdf[,3]>=PTHRESH[input$pv_thresh], 'cellRed', 'cellGreen'))
      #            )
      #)
     }
     tdf<-data.frame(tdf)
     tdf
  }, options=list(selection=list(target='row',selected=get_selection())))
    # HERE WE message OUT MICRORNA
    # messages table with selectable rows and conditional formatting
    output[[tbl2]] <- renderUI({
    tc<-NULL
     
 if (curtr==1)
   tbl<-input$testTbl1_1
 else
   tbl<-input$testTbl2_1
    selmir<-tbl[curtr]
    if (!is.null(tbl1))
      tc<-testcor(curtr)
    if (!is.null(tc)) {
      # restrict to top 10 mRNA
        tc<-na.omit(tc[1:10,])
        tc[,3]<-tc[,2]
        #tc[,4]<-tc[,2]
        colnames(tc)[3]<-colnames(tc)[2]
        colnames(tc)[2]<-"Link"

        mbnames<-tc
        mbnames[,1]<-gsub("_st","",mbnames[,1])
        mbnames[,1]<-gsub("[.]","_",mbnames[,1])
         mbnames[,1]<-gsub("_star","_5p",mbnames[,1])
        tc[,2]<-paste("<a href=http://www.mirbase.org/cgi-bin/query.pl?terms=",mbnames[,1]," target=_blank>","MiRBase","</a>",sep="")
         mbnames[,1]<-gsub("mir","miR",mbnames[,1])
         mbnames[,1]<-gsub("_","-",mbnames[,1])
	 #tc[,4]<-sapply(1:nrow(tc), function(x) {
	 	#cur_mm<-get.multimir(target=oldselmir1,mirna=mbnames[x,1])$validated$pubmedid;
	 #	cur_mm<-""
		#paste("<a href=www.ncbi.nlm.nih.gov/pubmed/",cur_mm[1],">",length(cur_mm),"</a>",sep="")
		 
	#})
	 #colnames(tc)[4]<-"Bind"
      HTML(df2html(tc, class = "tbl selRow", id = tbl2, 
                   cellClass = cbind(rep(NA, nrow(tc)), rep(NA, nrow(tc)), ifelse(abs(tc[,3])>=input$sp_thresh, 'cellGreen', 'cellRed'),NA)
         #          cellClass = cbind(rep(NA, nrow(tc)), rep(NA, nrow(tc)), ifelse(abs(tc[,3])>=input$sp_thresh, 'cellGreen', 'cellRed'))
                  )
      )
    }
  })
  
  curchoose_mir<-paste0("choose_mir",curtr)
output[[curchoose_mir]] <- renderUI({ 
     if (is.null(input$mirfile))
       return(NULL)

  selectInput("mir", "microRNA", colnames(mirInput()))
})

  # Compute the forumla text in a reactive expression since it is 
  formulaText1 <- reactive({
    selmir<-input$testTbl1_1
    if (selmir==null) {
     selmir = y[1,1]
    }
    paste(selmir," ~","pheno")
  })
  # Compute the forumla text in a reactive expression since it is 
  formulaText2 <- reactive({
    selmir<-input$testTbl2_1
    if (selmir==null) {
     selmir = y[1,1]
    }
    paste(selmir," ~","pheno")
  })

  # Generate a plot of the requested variable against x and only 
  # include outliers if requested
  curmirPlot<-paste0("mirPlot",curtr)
  output[[curmirPlot]] <- renderPlot({
    
    #curboxplot<-input$testTbl1_1_row_last_clicked
    par(mfrow=c(2,2))

     
    for (curboxplot in get_selection() ) {
    message("Doing boxplot")
    message(curboxplot)
    boxnames<-get_selection_names()
    if (curtr==1) {
      #message("DEBUG:pheno1")
      data=getPheno1()
      data=getPhenoExp(data,1,curboxplot)
  #FIXME concat selmir phenotype
      oldselmir<<-oldselmir1
    }
    else {
      #message("DEBUG:pheno2")
      data=getPheno2()
      data=getPhenoExp(data,2,curboxplot)
  #FIXME concat selmir phenotype
      oldselmir<<-oldselmir2
    }
    if (is.null(data)) {
 	return(NULL)
    }
    fold_thresh<-input$fold_thresh
    data$Group<-as.factor(data$Group)
   
    #fit = lm(Pheno~Group,data)
    y<-mirInput()
    y<-y[y[,1] %in% data[,1],]
    data<-data[data[,1] %in% y[,1],]
    y$Group<-data$Group
    y$Group<-as.factor(y$Group)
    nm<-names(y)[-ncol(y)]
    nm<-names(y)
    yaov<-0
    varlist<-nm
     #message("DEBUG:begin Anova step for boxplot")
    
    if (input$testtype=="Kruskal-Wallis")
    {  
       phenova<-kruskal.test(Pheno~Group,data)
       phenop<-phenova$"p.value"
    }
    else
    {  
      fit = lm(Pheno~Group,data)
      phenova<-anova(fit)
      phenop<-phenova$"Pr(>F)"[1]
    }
    whichbox<-which(boxnames$sels==curboxplot,)
    plotname<-boxnames[whichbox,]$names
    
    message(paste0("BOXPLOT DEBUG:data"))
    print(data)
    data$Group<-letters[as.numeric(as.character(data$Group))]
    
    boxplot(Pheno~Group,
	    data,
            outline = TRUE,
	    col = get_colors()[whichbox],
	    xlab = paste0(plotname," (p<",round(phenop,5),")"),
	    boxwex = 0.8
	  )
    abline(h=1,col='red')
    if (fold_thresh!=0) {
      abline(h=1+(fold_thresh/100),col='blue')
      abline(h=1-(fold_thresh/100),col='blue')
    }
     }},width=500,height=400) # this is the actual size of the plot output
#################

})
  } # done with both trees

  v <- reactiveValues(
    click1 = NULL,  # Represents the first mouse click, if any
    last = NULL, # last clicked/brushed plot
    range = NULL    # After two clicks, this stores the range of x
  )
   observeEvent(input$dblclick1_1, {
	v$click1=NULL
	v$range=NULL
	v$yrange=NULL
	v$last=1
  })
  observeEvent(input$dblclick1_2, {
	v$click1=NULL
	v$range=NULL
	v$yrange=NULL
	v$last=2
  })
  observeEvent(input$dblclick1_3, {
	v$click1=NULL
	v$range=NULL
	v$yrange=NULL
	v$last=3
  })
   observeEvent(input$plotbrush1_1, {
	v$range=range(input$plotbrush1_1$xmin,input$plotbrush1_1$xmax)
	v$yrange=range(input$plotbrush1_1$ymin,input$plotbrush1_1$ymax)
	v$last=1
  })
  observeEvent(input$plotbrush1_2, {
	v$range=range(input$plotbrush1_2$xmin,input$plotbrush1_2$xmax)
	v$yrange=range(input$plotbrush1_2$ymin,input$plotbrush1_2$ymax)
	v$last=2
  })
  observeEvent(input$plotbrush1_3, {
	v$range=range(input$plotbrush1_3$xmin,input$plotbrush1_3$xmax)
	v$yrange=range(input$plotbrush1_3$ymin,input$plotbrush1_3$ymax)
	v$last=3
  })

   observe ({
      #message(paste0("Table 5: ",  ifelse(is.null(input$testTbl1_2), "NULL", input$testTbl1_2 )))
      mrna_sel=input$testTbl1_2
      #message(input$testTbl1_2)
      message(numtests)
   })
  observe ({
       if (!is.null(input$testTbl1_1))
 	  oldselmir1<<-input$testTbl1_1
       if (!is.null(input$testTbl2_1))
 	  oldselmir2<<-input$testTbl2_1
       
     })

    # input$file1 will be NULL initially. After the user selects and uploads a 
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
    # columns. The 'datapath' column will contain the local filenames where the 
    # data can be found.

   output$downloadData <- downloadHandler(
   filename = function() {
     paste('data-', Sys.Date(), '.csv', sep='')
   },
   content = function(con) {
     write.csv(data, con)
   }
   )
   output$pcontents <- renderTable({
     datasetInput()
   }) 
   output$mircontents <- renderTable({
     mirInput()
   }) 
   output$mrnacontents <- renderTable({
     mrnaInput()
   }) 
    
messagesplit<-function(node,samples,varname,varval) {
  eqop="<"
  hdr=""  # should only be first split
  if (node%%2==1) { #right
    hdr=""
  }
  tmpstr<-paste0(hdr,"{\"samples\":",samples,",\"value\":[0],\"label\":\"","[",node,"]  ",varname,eqop,varval,"\",\"type\":\"split\",\"children\":[")
  return (tmpstr)
}

messageleaf<-function(node,samples,tr) {
  hdr=""
  clsr=""
  ds<-input[[paste0("depth_slider",tr)]]
  if ((node%%2)==0) { #right, because pheno labeling is 1-based
    hdr=""
    if (node>1) { #2,4,8
       clsrcnt<-0
       for (i in 1:ds) {
         clsrcnt<-clsrcnt+(node%%(2^i)==0) # add some for 2,4,8,16
        }
         for (i in 1:clsrcnt) {
           clsr<-paste0(clsr,"]}")
 	 }
       }
       else 
         clsr<-paste0(clsr,"]}")
   }
  #tmpstr<-paste0(hdr,"{\"samples\":",samples,",\"value\":[0],\"label\":\"","Pheno",node,"\"",",\"type\":\"leaf\"}",clsr)
  tmpstr<-paste0(hdr,"{\"samples\":",samples,",\"value\":[0],\"label\":\"","",letters[node],"\"",",\"type\":\"leaf\"}",clsr)
  return (tmpstr)
}

#   1
# 2   3
#45   67
#89abcdef

#124895ab36cd7ef
message_bst<-function(idx,sz) {
  tidx<-c(idx) 
  if (idx*2<sz) {
    tidx<-c(tidx,message_bst(idx*2,sz))
    tidx<-c(tidx,message_bst(idx*2+1,sz))
  }
  return(tidx)
}

find_best<-function() {
  
  if (input$findbest) {
    browser()
  }
  return(0)
}
tree_traverse<-function(x,tr) {
  
  retstr<-"" 
  #x <- datasetInput()
  #message("HERE IN tree_traverse")
  ds<-input[[paste0("depth_slider",tr)]]
  if (is.null(ds)) {
    return(retstr)
  }
  x<-x[c(colnames(x)[1],sapply(1:(2^ds-1),function(nx) input[[paste0("var",tr,"_",nx)]]))]
  if (is.null(x)) {
    return(retstr)
  }
  x <- x[complete.cases(x), ]
  curpheno<-x

  numnodes=2^(ds+1)
  tidx<-message_bst(1,numnodes)
  leaf_thresh=2^(ds)

  i<-1
  while (i <(numnodes-1)) {
    idx<-tidx[i]
    curpheno<-get_splitpheno(x,idx,tr,0,0)

    if (idx>=leaf_thresh) {
     phenonum<-(idx%%leaf_thresh)+1
     # left leaf
     retstr<-paste0(retstr,messageleaf(phenonum,nrow(curpheno),tr),",")
     i=i+1
     idx<-tidx[i]
     phenonum<-(idx%%leaf_thresh)+1
     curpheno<-get_splitpheno(x,idx,tr,0,0)
     retstr<-paste0(retstr,messageleaf(phenonum,nrow(curpheno),tr))
     # right leaf
     if (i<(numnodes-1)) {
     # up and over
       retstr<-paste0(retstr,",")
       i=i+1
       idx<-tidx[i]
     }
    }
    else {
      retstr<-paste(retstr,messagesplit(idx,nrow(curpheno),input[[paste0("var",tr,"_",idx)]],input[[paste0("range_slider",tr,"_",idx)]]),sep="")
      # traverse down
      i=i+1
      idx<-tidx[i]
    }
  }

  message (retstr)
  return(retstr[1])
}
########### end tree_traverse

for (tr in 1:NUMTREES) { # 2 trees
for (i in 1:n) {
    #message(paste0("DEBUG: In outer loop",tr," ",i))
local({
    curtr<-tr
    #make dynamic slider
    row <- i
    #message(paste0("DEBUG: In loop",curtr," ",row))
    curvar<-paste0("var",curtr,'_',row)
    #message("MAH in loopB4")
    curcol<-paste0('choose_columns',curtr,'_',row)
    output[[curcol]] <- renderUI({ 
       #message("DEBUG: render cur col")
       if (is.null(input$file1) && (input$loaddemo==0))
         {#message("MAH in null nputfile")
         return(NULL)
         }
	curvarname=''
	if (!is.null(input$deffile) || input$loaddemo>0)
        {
  	  curvarname<-subset(defInput(),Var==curvar)$Name
	  curvarname<-as.character(curvarname)
        }
    selectInput(curvar, paste0("Variable ",row),selected=curvarname,colnames(datasetInput()))
    })
    id<-paste0('range_slider',curtr,'_',row)
    output[[id]] <- renderUI({
    #message("DEBUG: render depth slider")
    depth<-input[[paste0("depth_slider",curtr)]]
      if (is.null(input$file1) && (input$loaddemo==0))
         return(NULL)
      if (is.null(input[[curvar]]))
  	 return(NULL)
    x <- datasetInput()
    x<-x[c(colnames(x)[1],sapply(1:(2^depth-1),function(nx) input[[paste0("var",curtr,"_",nx)]]))]
    #message(paste0("DEBUG:size of x is",nrow(x)))
    x <- x[complete.cases(x), ]
    #message(paste0("DEBUG:complete cases size of x is",nrow(x)))
    splitx<-get_splitpheno(x,row,curtr,0,0)
    if (!is.null(splitx) && (nrow(splitx)>0)) {
      splitx<-splitx[input[[curvar]]]
      slmin1 <- floor(min(splitx[input[[curvar]]]))
      slmax1 <- ceiling(max(splitx[input[[curvar]]]))
       name_idx=grep(input[[curvar]],colnames(splitx),fixed=TRUE)
       slmed1 <- ceiling(as.numeric(sapply(splitx[name_idx],median)))
       if (input$findbest) {
	  slmed1<-slmin1
	}
    }
    else {
      #message(paste0("DEBUG: Null when creating slider:",id))
      slmin1<-0
      slmax1<-0 
      slmed1<-0
    }
    if (!is.null(input$deffile) || input$loaddemo>0)
    {
  	  curvarsld<-subset(defInput(),Var==curvar)$Val
	  curvarsld<-as.numeric(as.character(curvarsld))
          slmed1<-curvarsld
    }
    else if (!is.null(input[[paste0("range_slider",curtr,"_",row)]]))
    {  
        	slmed1<-input[[paste0("range_slider",curtr,"_",row)]]
    }
    sliderInput(inputId = id,
                label = paste(""),
                #min = slmin1+(row<=(2^(depth-2)+1)), max = slmax1+(row>(2^(depth-2)+1)), value = slmed1)
                min = slmin1, max = slmax1+1, value = c(slmed1), animate=FALSE)
   })
   curhist<-paste0("hist",curtr,"_",row)
   curpline<-paste0("pline",curtr,"_",row)
   curdoofat<-paste0("doofat",curtr,"_",row)
   #slval<-input[[cursld]] 
   #input[[paste0("range_slider",curtr,"_",row)]]
   
   output[[curpline]] <- renderPlot({
	if (curtr==2) {} else {

   
   doallofat<-FALSE
   if (input[[curdoofat]] > olddoofat[curtr*row]) {
     olddoofat[curtr*row]<<-input[[curdoofat]]
     lastofat<-row
     ppvals<<-data.frame()
     ppsize<<-data.frame()
     doallofat<-TRUE
   }

   #output[[curpline]] <- renderChart2({
    #for (curboxplot in input$testTbl1_1_rows_selected) {
    curval<-input[[paste0("range_slider",curtr,"_",row)]]
     newx<-curval
    depth<-input[[paste0("depth_slider",curtr)]]
    x <- datasetInput()
    y <- mirInput()
    x<-x[c(colnames(x)[1],sapply(1:(2^depth-1),function(nx) input[[paste0("var",curtr,"_",nx)]]))]
    #message(paste0("DEBUG:size of x is",nrow(x)))
    x <- x[complete.cases(x), ]
    y<-y[y[,1] %in% x[,1],]
    x<-x[x[,1] %in% y[,1],]
    sels<-as.numeric(get_selection())+1 #add 1 because first column is id
    genes<-colnames(y)[sels]
    y<-y[genes] #add 1 because first column is patient id

    inputcolnames<-sub("X.","",colnames(data.frame(lapply(1:(2^depth-1),function(nx) input[[paste0("var1_",nx)]]))))
    inputcolnames<-sub("[.]","",inputcolnames)
    inputcurvar<-inputcolnames[row]   

    ofat<-data.frame()
    ds<-input$depth_slider1
    vals<-lapply(1:(2^ds-1),function(x) input[[paste0("range_slider",1,"_",x)]])
    message(vals)
    if (length(vals) < (2^ds-1))
      browser()
    if (doallofat) 
    {
      ofat<-get_ofat(x,y,depth,inputcurvar,curval,genes,doallofat)
      oldofat<<-ofat
    }
    else
      ofat<-oldofat

    minsize<-ofat$minsize
    ofat<-ofat[-ncol(ofat)] # strip off minsize info from ofat
    #message(paste0("DEBUG:complete cases size of x is",nrow(x)))
    splitx<-get_splitpheno(x,row,curtr,0,0)
      splitx<-splitx[inputcurvar]
      slmin1 <- floor(min(splitx[inputcurvar]))
      slmax1 <- ceiling(max(splitx[inputcurvar]))
    #}
     #ps<-get_pvals()
     par(mar=c(0,0.5,0,0))
     #par(pin=c(0.3,0.3))
     #par(new=T)
     
     #newx<-t(t(rep(newx,length(newy))))
     newx<-as.numeric(rownames(ofat))
     newx<-rep(newx,ncol(ofat))
     newc<-get_colors()[col(ofat)]
     newy<-unlist(ofat)
     names(newy)=c()
     newv<-inputcurvar
     pindex<-row
     newd<-unlist(lapply(minsize,function(x) if(x>4) 0.9 else (2*x)/10))
     newpoints<-data.frame(pindex,newv,newx,newy,newc,newd)
     message(paste0("curval:",curval))
     message(paste0("newpoints:",newpoints))
     message(paste0("ppvals:",ppvals))
     message(paste0("ppsize:",ppsize))
     message(paste0("newv:",newv))
     message(paste0("minsize:",minsize))
     ppvals<<-rbind(ppvals,newpoints)
     #ppvals<<-cbind(ppvals,minsize)
     curpvals<-ppvals[ppvals$pindex==pindex,]
     cxpoints<-curpvals[curpvals$newv==newv,]$newx
     cypoints<-curpvals[curpvals$newv==newv,]$newy
     cxcols<-curpvals[curpvals$newv==newv,]$newc
     dotsizes<-curpvals[curpvals$newv==newv,]$newd
     cxcols<-adjustcolor(cxcols,alpha.f=0.5)
     plxlim<-v$range
     plylim<-v$yrange
     if (is.null(plxlim) || row != v$last) {
     plxlim<-c(slmin1,slmax1)
     }
     if (is.null(plylim) || row != v$last) {
     plylim<-c(0,6)
     }
     message(paste0("dotsizes:",dotsizes))
     plot(cxpoints,-log10(cypoints), axes=FALSE,xlab='',ylab='',col = cxcols,ps=2,cex=dotsizes,cex.axis=0.6,pch=19,xlim=plxlim,ylim=plylim)
     axis(side=1,cex.axis=0.6)
     if (input[[paste0("ptrend",curtr,"_",row)]]) {
	linedf<-data.frame(cxpoints,-log10(cypoints),cxcols)
	linedf<-linedf[order(cxpoints),]
	for (lcol in unique(cxcols)) {
          lines(linedf[linedf$cxcols==lcol,][1:2],col=lcol)
	}
       #lines(loess(-log10(cypoints)~cxpoints,degree=2,span=0.5),col=get_colors())
    }   
     
     #gplotd<-data.frame(xv=cxpoints,yv=-log10(cypoints))
     #p<-ggplot(gplotd,aes(x=xv,y=yv))+geom_point(aes(size=10))
     #p<-p+theme(axis.line=element_blank(),axis.text.x=element_blank(),
     #     axis.text.y=element_blank(),axis.ticks=element_blank(),
     #     axis.title.x=element_blank(),
     #     axis.title.y=element_blank())
     #p<-p+ylim(0,4) 
     #p<-p+xlim(slmin1,slmax1) 
     print(p)
     abline(v=curval)
     abline(h=c(1,2,3),col="gray")
     #points(newx,newy)
    # lines(ps)
   
   }},width=200,height=50)
   output[[curhist]] <- renderPlot({
     #message("DEBUG: render hist Plot")
     x <- datasetInput()
     depth<-input[[paste0("depth_slider",curtr)]]
     x<-x[c(colnames(x)[1],sapply(1:(2^depth-1),function(nx) input[[paste0("var",curtr,"_",nx)]]))]
     x <- x[complete.cases(x), ]
     splitx<-get_splitpheno(x,row,curtr,0,0)
     curnode<-splitx
     if (is.null(curnode) || nrow(curnode)==0) {
        #message("DEBUG: Null histogram")
     } else {
       bp<-curnode[input[[curvar]]]
     if (length(bp)>0) {
       par(pin=c(0.3,0.3))
       par(mar=c(0,0.5,0,0))
       normals[curhist]<<-shapiro.test(unlist(bp))$p.value
       hist(bp[,1],main=NULL,col="blue",xlab=NULL,ylab=NULL,axes=FALSE,right=FALSE)
     }
    }
   },width=200,height=20)
})
}
}


observe({
   lvls<-input$depth_slider1
   message(paste0("DEBUG:depth_slider observe:",length(oldpsy)))
})

observe({

if (length(input$testTbl1_1_rows_selected)>0) {
	selrow1<<-as.numeric(input$testTbl1_1_rows_selected)
	message(paste0("DEBUG:Selected rows observe:",selrow1))
}
})

observe({
#message("DEBUG:range_slider observe")

for (t in 1:2) {
for (i in 1:n) {
  input[[paste0("range_slider",t,"_",i)]]
}}
for (i in 1:n) {
  #message(input[[paste0("range_slider1_",i)]])
}
input$range_slider1_1 # Do take a dependency on input$saveButton
input$range_slider1_2 # Do take a dependency on input$saveButton
input$range_slider1_3 # Do take a dependency on input$saveButton


# isolate a whole block
data <- isolate({
#data <- ({
lvls<-input$depth_slider1
lvls2<-input$depth_slider2
#message("IN data block")
a <- input$range_slider1_1
b <- input$range_slider1_2
c <- input$range_slider1_3
d <- input$range_slider2_1
e <- input$range_slider2_2
f <- input$range_slider2_3
x <- datasetInput()


pattern<-tree_traverse(x,1)
sink("./www/tmp.json")
cat(pattern)
sink()

file.rename("./www/tmp.json",paste0("./www/tree",thetime,"_1.json"))

pattern<-tree_traverse(x,2)
sink("./www/tmp.json")
cat(pattern)
sink()
file.rename("./www/tmp.json",paste0("./www/tree",thetime,"_2.json"))

})
})
})
