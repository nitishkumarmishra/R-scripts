library(GDCRNATools)
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/', "/")
mirdir <- paste(project, 'miRNAs', sep='/', "/")

dir.create(project)
dir.create(rnadir)
dir.create(mirdir)
setwd(rnadir)

gdcRNADownload(project.id = 'TCGA-CHOL', data.type = 'RNAseq', write.manifest = TRUE, directory = rnadir)
#gdcRNADownload(project.id     = 'TCGA-CHOL', data.type      = 'miRNAs',  write.manifest = TRUE, directory  = "../miRNAs/")

# Manually use manifest file name in next line
gdcRNADownload(manifest = "TCGA-CHOL.RNAseq.gdc_manifest.2018-07-26T10-23-05.txt")

metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',data.type  = 'RNAseq',write.meta = FALSE)

rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, data.type = 'RNAseq', path="Data", organized = FALSE)
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR, data.type = 'miRNAs', path="Data", organized = FALSE)
#mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,  path      = mirdir, data.type = 'miRNAs')


setwd("../miRNAs/")
gdcRNADownload(project.id  = 'TCGA-CHOL', data.type = 'miRNAs',  write.manifest = TRUE)
grepl("TCGA", list.files())

#gdcRNADownload(manifest = grepl("TCGA", list.files()))
gdcRNADownload(manifest = "TCGA-CHOL.miRNAs.gdc_manifest.2018-07-26T11-01-52.txt")
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR, data.type = 'miRNAs', path="Data", organized = FALSE)


rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = TRUE)
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
#mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = TRUE)
rnaCounts <- rnaCounts[rownames(rnaExpr),] ## Select after filter
mirCounts <- mirCounts[rownames(mirExpr),]##

DEGAll <- gdcDEAnalysis(counts = rnaCounts, group = metaMatrix.RNA$sample_type,comparison = 'PrimaryTumor-SolidTissueNormal', method = 'limma')
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

### ceRNAs network analysis using internal databases
## I didn't find network if I am using filter=TRUE in gdcVoomNormalization
ceOutput <- gdcCEAnalysis(lnc = rownames(deLNC), pc = rownames(dePC),lnc.targets = 'starBase', pc.targets  = 'starBase', rna.expr = rnaExpr, mir.expr = mirExpr)

### ceRNAs network analysis using user-provided datasets
# load miRNA-lncRNA interactions
data(lncTarget)
data(pcTarget)
ceOutput <- gdcCEAnalysis(lnc = rownames(deLNC), pc = rownames(dePC), lnc.targets = lncTarget, pc.targets  = pcTarget, rna.expr = rnaExpr, mir.expr = mirExpr)

### Network visulization in Cytoscape

# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]


# Edges and nodes can be simply imported into Cytoscape 
# for network visualization
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges') edges[1:5,]

nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
write.table(edges, file='edges.txt', sep='\t', quote=F) ### Network of Cytoscape
write.table(nodes, file='nodes.txt', sep='\t', quote=F) ### Table of Cytoscape


######################
### Correlation plot on a local webpage

shinyCorPlot(gene1  = rownames(deLNC), gene2  = rownames(dePC), rna.expr = rnaExpr, metadata = metaMatrix.RNA)

# CoxPH analysis
survOutput <- gdcSurvivalAnalysis(gene = rownames(deALL), method   = 'coxph', rna.expr = rnaExpr, metadata = metaMatrix.RNA)


# KM analysis
survOutput <- gdcSurvivalAnalysis(gene = rownames(deALL), method = 'KM', rna.expr = rnaExpr, metadata = metaMatrix.RNA, sep = 'median')

# KM plot on a local webpage by shinyKMPlot
shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, metadata = metaMatrix.RNA)

## Enrichment analysis
enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE, level = 3)## level3 will be fast

# Barplot
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

# Bubble plot
gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)


# View pathway maps on a local webpage
library(pathview)

deg <- deALL$logFC
names(deg) <- rownames(deALL)
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])

shinyPathview(deg, pathways = pathways, directory = 'pathview')
