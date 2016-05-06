library(biomaRt)
library(GOstats)
library(topGO)


listMarts()
ensembl = useMart("ensembl", dataset="scerevisiae_gene_ensembl")
data<-listDatasets(ensembl)
filt<-listFilters(ensembl)
attr<-listAttributes(ensembl)
query="CDC19"
testlst<-c("CDC15", "CDC16", "CDC17", "CDC18", "CDC19") 
filter<-"go_id"
filter2<-"wikigene_name"
b<-getBM(attributes=c("entrezgene", "go_id", "name_1006"), filters=filter2, values=lst, mart=ensembl)
e<-getBM(attributes=c("entrezgene", "description"), filters=filter2, values=query, mart=ensembl)
testgo<-getBM(attributes=c(filter2, "go_id"), filters=filter2, values=testlst, mart=ensembl)
geneToGo<-getBM(attributes=c("wikigene_name", "go_id"), filters="with_go_go", values=TRUE, mart=ensembl)

####Format Gene to Go list!
gene.to.GO <- split(geneToGo$go_id, geneToGo$wikigene_name)
geneToGO <- lapply(gene.to.GO, unique)  # to remove duplicates

write.csv(geneToGO, file = "sacc_list.csv", na="NA")

#####################

lst<-read.csv("/Users/nicholaswaters/Documents/GU_R/biomaRT/testGoGene.csv")
go<-lst$GO.ID
pval<-lst$p.value
gen<-lst$List.of.Genes
filter2<-"wikigene_name"
lizGenes<-

g<-getBM(attributes="entrezgene", filters=filter2, values=gen, mart=ensembl)

##   Look it F$!@ing works!!!!!!
GOdata <- new("topGOdata", ontology = "MF", allGenes = testList, annot = annFUN.gene2GO, gene2GO = geneToGO)


##but this doesnt work
tgd<-new("topGOdata", ontology="CC", allGenes=geneList, nodeSize=1, description = "Ensembl GO enrichment", 
         annot = annFUN.gene2GO, gene2GO =testgo)
######################

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

