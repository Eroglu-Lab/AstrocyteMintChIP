
library(Gviz)
library(GenomicRanges)
library(Mus.musculus)
library(rtracklayer)
library(stringr)
library(shiny)
library(glue)
library(coMET)
library(shinythemes)
library(biomaRt)
library(ggplot2)
library(viridis)
library(aws.s3)

#read in data#
Sys.setenv("AWS_ACCESS_KEY_ID" = "AKIA2KPURUMOTSR5SJ6N",
           "AWS_SECRET_ACCESS_KEY" = "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",
           "AWS_DEFAULT_REGION" = "us-east-2")

#read in data#
Bucket <- get_bucket('mintchip-dev-bigwigs')
objectList <- c()
for (i in 1:length(Bucket)){
  objectList[i] <- Bucket[[i]]$Key
}

bwList <- list()
for (i in 1:length(objectList)){
  bwList[i] <- s3read_using(object=objectList[i], bucket='mintchip-dev-bigwigs', FUN=function(x) import.bw(x, as='GRanges'))
}
names(bwList) <- gsub('.bam.bw', '', objectList)
names(bwList) <- gsub('.bw', '', names(bwList))


refGene <- read.csv('mm10_refGene_noDup.csv')

fpkm <- read.csv('fpkmLongFormat.csv')
fpkm$age <- factor(fpkm$age, levels=c('P7', 'P14', 'P21', 'P150'))

### define function to obtain genome coords for input gene ###
getCoordinates <- function(File, genename){
  x <- File[File$geneName == genename, ]
  chr <- x[,1]
  start <- x[,4]
  end <- x[,5]
  chr2 <- as.numeric(str_match(chr, 
                               '(\\d{1,})')[,2])
  coords <- list(chr, start, end, chr2)
}

hmarks <- c('H3K4me1', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K27me3', 'H3K36me3')
#ages <- c('P7', 'P14', 'P21', 'P150')
