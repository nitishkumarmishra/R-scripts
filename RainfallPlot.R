library(naturalsort);library(gtrellis);library(circlize);library(ComplexHeatmap)
#load("DiffMeth.RData")
### Results is limma differential methylation results file from DiffMeth.RData
temp <- results[,c("CHR", "MAPINFO", "MAPINFO", "log2FC")]
temp$CHR <- paste("chr", temp$CHR, sep = "")
### change colname, add start, end, foldchange in dataframe
names(temp)[1] <- c("chr"); names(temp)[2] <- c("start"); names(temp)[3] <- c("end")
temp$direction <- ifelse(temp$log2FC >0, "Hyper", "Hypo")
temp$log2FC <- NULL; df <- temp
df <- df[naturalorder(df$chr),]### Order dataframe with chromosome name
DMR_hyper = df[df$direction == "Hyper", ]
DMR_hypo = df[df$direction == "Hypo", ]
DMR_hyper_density = genomicDensity(DMR_hyper, window.size = 2000000)
DMR_hypo_density = genomicDensity(DMR_hypo, window.size = 2000000)
max_density = max(c(DMR_hyper_density[[4]], DMR_hypo_density[[4]]))
DMR_hyper_rainfall = rainfallTransform(DMR_hyper)
DMR_hypo_rainfall = rainfallTransform(DMR_hypo)
cm = ColorMapping(levels = c("Hyper", "Hypo"), colors = c("red", "blue4"))
lgd = color_mapping_legend(cm, title = "Direction", plot = FALSE)
### I have used all differentially methylated CpGs with p < 0.01
pdf("RainfallPlot.pdf", width = 10, height = 10)
#pdf("RainfallPlotDeltabeta-0.1.pdf", width = 10, height = 10)#with delata beta 0.1
#gtrellis_layout(category = c(paste0("chr", 1:22), c("chrX", "chrY")), n_track = 3, ncol = 4, byrow = FALSE,
gtrellis_layout(category = paste0("chr", 1:22), n_track = 3, ncol = 4, byrow = FALSE,
                track_axis = TRUE, 
                track_height = c(1, 0.5, 0.5), 
                track_ylim = c(0, 8, 0, max_density, 0, max_density),
                track_ylab = c("log10(dist)", "Hyper", "Hypo"),
                add_name_track = TRUE, add_ideogram_track = TRUE,
                legend = lgd, title = "Hyper and hypomethylated DM CpGs")
### plot for Hyper and Hypo CpGs with log10(dist)
add_points_track(DMR_hyper_rainfall, log10(DMR_hyper_rainfall[[4]]),
                 pch = 16, size = unit(0.9, "mm"), gp = gpar(col = "red"))
add_points_track(DMR_hypo_rainfall, log10(DMR_hypo_rainfall[[4]]), track = current_track(),
                 pch = 16, size = unit(0.9, "mm"), gp = gpar(col = "blue4"))
# track for genomic density
add_lines_track(DMR_hyper_density, DMR_hyper_density[[4]], area = TRUE, 
                gp = gpar(fill = "red", col = NA))
add_lines_track(DMR_hypo_density, DMR_hypo_density[[4]], area = TRUE,
                gp = gpar(fill = "blue4", col = NA))
dev.off()


rm( commonID, Meth, GeneExp, contrast.matrix, Clinic, methFit, methFit2, methRresults, methRresults1, Tumor.impue, Normal.impue, Tumor.impue.BMIQ, Tumor.Meth, Normal.impue.BMIQ, 
    Normal.Meth, c1, c2, t.tumor, probe.features, probeInfoALL.lv, design, BMIQ.Meth, BMIQ.Meth.M, myNorm, numNAs, limma.results)
save.image("rainfall.RData")
save.image("rainfallDelat0.1.RData")#with delata beta 0.1