server <- function(input, output) {

### receive information about the input gene ###
chr <- reactive({ getCoordinates(refGene, input$Gene)[[1]] })
Start <- reactive({ getCoordinates(refGene, input$Gene)[[2]] })
End <- reactive({ getCoordinates(refGene, input$Gene)[[3]] })
chrNum <- reactive({ getCoordinates(refGene, input$Gene)[[4]] })

### make annotation tracks ### 
genesUcsctrack <- reactive({ 
  bgrTrack <- BiomartGeneRegionTrack(genome="mm10",
                                     start=Start(),
                                     end=End(),
                                     chromosome = chrNum(),
                                     biomart=useEnsembl(biomart = 'genes', 
                                                        dataset = 'mmusculus_gene_ensembl',
                                                        version = 101),
                                     name="ENSEMBL",
                                     fontcolor.title='black',
                                     fontsize=16,
                                     fontsize.group=20) })

  axisTrack <- GenomeAxisTrack()

### make data tracks ###

k27m_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K27me3\nP7', col.histogram="grey70", fill.histogram="grey70", 
            fontsize=16, fontcolor.title='black', background.title='grey70', ylim=c(0,100)) })
k27m_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K27me3\nP14', col.histogram="grey40", fill.histogram="grey40", 
            fontsize=16, fontcolor.title='black', background.title='grey40',ylim=c(0,100)) })
k27m_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K27me3\nP21', col.histogram="grey20", fill.histogram="grey20", 
            fontsize=16, fontcolor.title='white', background.title='grey20',ylim=c(0,100)) })
k27m_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K27me3\nAdult', col.histogram="#000000", fill.histogram="#000000", 
            fontsize=16, fontcolor.title='white', background.title='#000000',ylim=c(0,100)) })

k27_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K27ac\nP7', col.histogram="chocolate1", fill.histogram="chocolate1", 
            fontsize=16, fontcolor.title='black', background.title='chocolate1',ylim=c(0,200)) })
k27_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K27ac\nP14', col.histogram="chocolate2", fill.histogram="chocolate2", 
            fontsize=16, fontcolor.title='black', background.title='chocolate2',ylim=c(0,200)) })
k27_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K27ac\nP21', col.histogram="chocolate3", fill.histogram="chocolate3", 
            fontsize=16, fontcolor.title='black', background.title='chocolate3',ylim=c(0,200)) })
k27_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K27ac\nAdult', col.histogram="chocolate4", fill.histogram="chocolate4", 
            fontsize=16, fontcolor.title='black', background.title='chocolate4',ylim=c(0,200)) })

k36_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K36me3\nP7', col.histogram="lightpink1", fill.histogram="lightpink1", 
            fontsize=16, fontcolor.title='black', background.title='lightpink1',ylim=c(0,200)) })
k36_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K36me3\nP14', col.histogram="lightpink2", fill.histogram="lightpink2", 
            fontsize=16, fontcolor.title='black', background.title='lightpink2',ylim=c(0,200)) })
k36_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K36me3\nP21', col.histogram="lightpink3", fill.histogram="lightpink3", 
            fontsize=16, fontcolor.title='black', background.title='lightpink3',ylim=c(0,200)) })
k36_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K36me3\nAdult', col.histogram="lightpink4", fill.histogram="lightpink4", 
            fontsize=16, fontcolor.title='black', background.title='lightpink4', ylim=c(0,200)) })

k9_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K9me3\nP7', col.histogram="steelblue1", fill.histogram="steelblue1", 
            fontsize=16, fontcolor.title='black', background.title='steelblue1', ylim=c(0,100)) })
k9_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K9me3\nP14', col.histogram="steelblue2", fill.histogram="steelblue2", 
            fontsize=16, fontcolor.title='black', background.title='steelblue2',ylim=c(0,100)) })
k9_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K9me3\nP21', col.histogram="steelblue3", fill.histogram="steelblue3", 
            fontsize=16, fontcolor.title='black', background.title='steelblue3',ylim=c(0,100)) })
k9_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K9me3\nAdult', col.histogram="steelblue4", fill.histogram="steelblue4", 
            fontsize=16, fontcolor.title='black', background.title='steelblue4',ylim=c(0,100)) })

k41_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K4me1\nP7', col.histogram="orangered1", fill.histogram="orangered1", 
            fontsize=16, fontcolor.title='black', background.title='orangered1',ylim=c(0,200)) })
k41_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K4me1\nP14', col.histogram="orangered2", fill.histogram="orangered2", 
            fontsize=16, fontcolor.title='black',background.title='orangered2', ylim=c(0,200)) })
k41_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K4me1\nP21', col.histogram="orangered3", fill.histogram="orangered3", 
            fontsize=16, fontcolor.title='black', background.title='orangered3',ylim=c(0,200)) })
k41_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K4me1\nAdult', col.histogram="orangered4", fill.histogram="orangered4", 
            fontsize=16, fontcolor.title='black', background.title='orangered4',ylim=c(0,200)) })

k43_7 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P7')]], type='histogram',
            name='H3K4me3\nP7', col.histogram="aquamarine1", fill.histogram="aquamarine1", 
            fontsize=16, fontcolor.title='black', background.title='aquamarine1',ylim=c(0,200)) })
k43_14 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P14')]], type='histogram',
            name='H3K4me3\nP14', col.histogram="aquamarine2", fill.histogram="aquamarine2", 
            fontsize=16, fontcolor.title='black', background.title='aquamarine2',ylim=c(0,200)) })
k43_21 <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'P21')]], type='histogram',
            name='H3K4me3\nP21', col.histogram="aquamarine3", fill.histogram="aquamarine3", 
            fontsize=16, fontcolor.title='black',background.title='aquamarine3', ylim=c(0,200)) })
k43_a <- reactive({ 
  DataTrack(range=bwList[[paste0(input$mark, '_', 'Adult')]], type='histogram',
            name='H3K4me3\nAdult', col.histogram="aquamarine4", fill.histogram="aquamarine4", 
            fontsize=16, fontcolor.title='black', background.title='aquamarine4',ylim=c(0,200)) })

## chromosome track ##
itrack <- reactive ({ IdeogramTrack(genome = 'mm10', chromosome = chr(), fontsize=20) })

dataTrackLists <- reactive ({ 
  list(k27m_7(),k27m_14(),k27m_21(),k27m_a(),
       k27_7(),k27_14(),k27_21(),k27_a(),
       k36_7(),k36_14(),k36_21(),k36_a(),
       k9_7(),k9_14(),k9_21(),k9_a(),
       k41_7(),k41_14(),k41_21(),k41_a(),
       k43_7(),k43_14(),k43_21(),k43_a()) })

output$Summary <- renderText({ 
                 'Welcome to the browser of developmental astrocyte histone modifications.
                 Data were collected from P7, P14, P21 and P150 (Adult) astrocyte cortices
                 In general the histone modifications map to the following genomic features:
                 H3K4me1: Active enhancers
                 H3K4me3: Active promoters / TSS
                 H3K9me3: Repressed regions
                 H3K27ac: Active promoters / enhancers
                 H3K27me3: Repressed gene body
                 H3K36me3: Active transcription, gene body'
                
                 })

output$Tracks <- 
  renderPlot({ 
  if ('H3K4me1' %in% input$mark){

    
    plotTracks(list(k41_7(), k41_14(), k41_21(), k41_a(),genesUcsctrack(),itrack(), axisTrack), 
           chromosome=chrNum(),
           from=Start() - 5000, 
           to=End() + 5000,
           transcriptAnnotation="symbol",
           geneSymbol=T)

  } else if ('H3K4me3' %in% input$mark){
  
    plotTracks(list(k43_7(), k43_14(), k43_21(), k43_a(),genesUcsctrack(),itrack(), axisTrack), 
    chromosome=chrNum(),
    from=Start() - 5000, 
    to=End() + 5000,
    transcriptAnnotation="symbol",
    geneSymbol=T)
    
  } else if ('H3K9me3' %in% input$mark){
    
    plotTracks(list(k9_7(), k9_14(), k9_21(), k9_a(),genesUcsctrack(),itrack(), axisTrack), 
               chromosome=chrNum(),
               from=Start() - 5000, 
               to=End() + 5000,
               transcriptAnnotation="symbol",
               geneSymbol=T)
    
  } else if ('H3K36me3' %in% input$mark){
    
    plotTracks(list(k36_7(), k36_14(), k36_21(), k36_a(),genesUcsctrack(),itrack(), axisTrack), 
               chromosome=chrNum(),
               from=Start() - 5000, 
               to=End() + 5000,
               transcriptAnnotation="symbol",
               geneSymbol=T)
    
  } else if ('H3K27me3' %in% input$mark){
    
    plotTracks(list(k27m_7(), k27m_14(), k27m_21(), k27m_a(),genesUcsctrack(),itrack(), axisTrack), 
               chromosome=chrNum(),
               from=Start() - 5000, 
               to=End() + 5000,
               transcriptAnnotation="symbol",
               geneSymbol=T)
    
  } else if ('H3K27ac' %in% input$mark){
    
    plotTracks(list(k27_7(), k27_14(), k27_21(), k27_a(),genesUcsctrack(),itrack(), axisTrack), 
               chromosome=chrNum(),
               from=Start() - 5000, 
               to=End() + 5000,
               transcriptAnnotation="symbol",
               geneSymbol=T)
}

})

######FPKM plots
fpkm_r <- reactive({
  fpkm %>%
    dplyr::filter(genes == input$Gene)})

output$rna <- renderPlot({
  
  ggplot(fpkm_r(), aes(x=age, y=fpkm.long, fill=age)) + geom_bar(stat='identity') + 
    theme_minimal() + ylab('FPKM') +  xlab('') + scale_fill_manual(values=c('gray75', 'gray50', 'gray25', 'gray0')) + ggtitle(print(input$gene))

})


}
