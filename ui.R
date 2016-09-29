library(shiny)
library(shinydashboard)
library(parallel)
library(miscTools)
#library(multiMiR)
library(rCharts)
library(DT)
source('df2html.R')
tr=1

tree1_string=paste0("<html><body><script type=\"text/javascript\">var treename = \"./tree",thetime,"_1.json\"\nvar tree1 = PhenoTree(treename)\n</script></body></html>")
tree2_string=paste0("<html><script type=\"text/javascript\">var tree2name = \"./tree",thetime,"_2.json\"\nvar tree2 = PhenoTree(tree2name)\n</script></html>")

shinyUI(fluidPage(
  titlePanel("PEAX"),
  tags$head(
     tags$style(HTML("
       .skin-blue .main-sidebar {
       		background-color: #ffffff;
                     .sidebar { height: 90vh; overflow-y: scroll; }
                     },
                     " )
     )
   ),
tags$style(type="text/css",
 ".shiny-output-error { visibility: hidden; }",
".shiny-output-error:before { visibility: hidden; }"
),
  sidebarLayout(
  #sidebarPanel(style = "background-color: #ffffff;width:250px", 
  sidebarPanel(style = "background-color: #ffffff;",width=2, 
tabsetPanel(id="leftTabs", style = "background-color: #ffffff;overflow-y:scroll;",


tabPanel(style = "background-color: #ffffff;", "Data",
actionButton("loaddemo","Load Demo Data"),
  fileInput('file1', 'Phenotype File',
            accept=c('text/csv', 'text/comma-separated-values,text/plain')),
  fileInput('mirfile', 'Expression File 1',
            accept=c('text/csv', 'text/comma-separated-values,text/plain')),
  fileInput('mrnafile', 'Expression File 2',
            accept=c('text/csv', 'text/comma-separated-values,text/plain')),
  downloadLink('downloadFile', 'Download Tree'),
  fileInput('evidence', 'miR-mRNA Evidence',
            accept=c('text/csv', 'text/comma-separated-values,text/plain')),
  checkboxInput('header', 'Header', TRUE),
  checkboxInput('transpose', 'Transpose', FALSE),
  radioButtons('sep', 'Separator',
               c(Comma=',',
                 Semicolon=';',
                 Tab='\t'),
               'Comma'),
  radioButtons('quote', 'Quote',
               c(None='',
                 'Double Quote'='"',
                 'Single Quote'="'"),
               'Double Quote'
 )),
tabPanel(style = "background-color: #ffffff;", id="Analysis","Analysis",
mclapply(1:n, function(i) {
conditionalPanel(condition=paste0("input.depth_slider",tr,">",log2(i)),
   div(class="row-fluid",htmlOutput(paste0("choose_columns",tr,"_",i))),
   div(class="row-fluid",plotOutput(paste0("hist",tr,"_",i)), style="height:20px"),
   div(style="display:inline-block; vertical-align: middle;",class="row-fluid",checkboxInput(paste0('ptrend',tr,"_",i),"Trend", FALSE)),
   div(style="display:inline-block; vertical-align: middle;",actionButton(paste0('doofat',tr,"_",i),"OFAT")),
   div(class="row-fluid",plotOutput(paste0("pline",tr,"_",i),clickId=paste0("plotclick",tr,"_",i),dblclick=paste0("dblclick",tr,"_",i),
		brush=paste0("plotbrush",tr,"_",i)), style="height:50px"),
   div(class="row-fluid",uiOutput(paste0("range_slider",tr,"_",i)))
    )})
)

)

),
mainPanel(
#includeCSS("www/custom.css"),
  conditionalPanel(condition="input.leftTabs=='Data'",
    tabsetPanel(id="tabSelected",
  	tabPanel(style = "background-color: #ffffff;", "Phenotype Data",tableOutput('pcontents')),
  	tabPanel(style = "background-color: #ffffff;", "Exp1 Data",tableOutput('mircontents')),
  	tabPanel(style = "background-color: #ffffff;", "Exp2 Data",tableOutput('mrnacontents')))),
conditionalPanel(condition="input.leftTabs=='Analysis'",
#column(12,
   wellPanel(style = "background-color: #ffffff;",
   div(style="display:inline-block; vertical-align: middle;",fileInput('deffile', 'Load Tree',
            accept=c('text/csv', 'text/comma-separated-values,text/plain'))),
   div(style="display:inline-block; vertical-align: middle;",sliderInput("depth_slider1", "Tree Depth", min=1, max=4, value=2)),
   div(style="display:inline-block; vertical-align: middle;", sliderInput("pv_thresh","p-value threshold",min=1,max=7,value=1,step=1,ticks=TRUE)),    
   div(style="display:inline-block; vertical-align: middle;", sliderInput("sp_thresh","Correlation threshold",min=0.4,max=1,value=0.6,step=0.1)),
   div(style="display:inline-block; vertical-align: middle;", sliderInput("fold_thresh","Fold-change threshold",min=0,max=100,value=0,step=10)),

   div(style="display:inline-block; vertical-align: middle; horizontal-align: middle",
   
   textOutput("normality"),
  radioButtons('testtype', '                           ',
               c('Kruskal-Wallis',
                 'Anova'),
               inline=FALSE,
               'Kruskal-Wallis'),

   textOutput("numtests"))),     
      fluidRow(
	 column(4,box(includeHTML("www/tree.html"),width=NULL,
            HTML(tree1_string),height="600px")),
       column(3,
	
      radioButtons('selector', 'Selection',
               c('Table'='table',
                 'Volcano'='volcano'),
		inline=TRUE),
      conditionalPanel(condition="input.selector=='table'",
      	actionButton("search","Find Trend"),
      	actionButton("searchall","Find All"),
      	actionButton("findbest","Find Best"),
      	DT::dataTableOutput(outputId="testTbl1_1"),uiOutput(outputId="testTbl1_2"),width=NULL),
      conditionalPanel(condition="input.selector=='volcano'")),
      column(5,
            box(plotOutput("mirPlot1",height="600px",width="600px"),width=NULL)))

)))
)) # shinyUI,fluidPage,fluidRow
