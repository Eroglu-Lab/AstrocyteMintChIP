#start app#
ui <- fluidPage(theme=shinytheme('united'),
  titlePanel("Astrocyte developmental histone modifications browser"),
  
  sidebarLayout(
    sidebarPanel(selectInput("Gene", 
                             label = "Choose a gene",
                             choices = refGene$geneName,
                             selected = "Gfap"),
                 radioButtons("mark", "Select histone mark(s)", choices=hmarks, selected='H3K36me3' )),
    
    
    mainPanel(
              tabsetPanel(type = "tab",
                          tabPanel("Information", verbatimTextOutput("Summary")),
                          tabPanel("Histone Modification Tracks", plotOutput('Tracks')),
                          tabPanel("RNAseq FPKM", plotOutput('rna'))
                          )
  ))
)
