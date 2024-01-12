############################################################################
######### DE and WGCNA after accounting for cell-type abundance ############
############################################################################


# Install all required packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("lumiHumanAll.db")

#For DE
library(lumi)
library(lumiHumanIDMapping) 
library(lumiHumanAll.db) 
library(sva) 
library(limma)
library(annotate)
library(biomaRt)
library(siggenes)
library(broom) 
library(WGCNA)
library(tools)
library(ggplot2)

#String operations
library(stringr)

#Reading and writing tables
library(readr)
library(openxlsx)

#Plotting
library(ggplot2)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Data arrangement
library(dplyr)
library(tidyr)

#Functional programming
library(magrittr)
library(purrr)
library(rlist)

#setworkingdirectory
#setwd("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/regressed_pipeline/")
setwd("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/regressed_pipeline/rm_batch5_subs")


### FUNCTIONS ###
#Histogram
Histogram <- function(filename, lumi.object) {
  expr.df <- exprs(lumi.object) %>% t %>% data.frame
  dataset.addvars <- mutate(expr.df, Sample.Name = sampleNames(lumi.object), Status = lumi.object$Diagnosis22q)
  dataset.m <- gather(dataset.addvars, nuID, Expression, -Sample.Name, -Status)
  
  p <- ggplot(dataset.m, aes(Expression, group = Sample.Name, col = factor(Status))) + geom_density() + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
  CairoPDF(filename, height = 5, width = 9)
  print(p)
  dev.off()
}

#Samplewise connectivity plot
ConnectivityPlot <- function(filename, dataset, maintitle) {
  norm.adj <- (0.5 + 0.5 * bicor(exprs(dataset)))
  colnames(norm.adj) <- dataset$Slide
  rownames(norm.adj) <- dataset$Slide
  net.summary <- fundamentalNetworkConcepts(norm.adj)
  net.connectivity <- net.summary$Connectivity
  connectivity.zscore <- (net.connectivity - mean(net.connectivity)) / sqrt(var(net.connectivity))
  
  connectivity.plot <- data.frame(Slide.ID = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
  p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Slide.ID) )
  p <- p + geom_text(size = 4, colour = "red")
  p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
  CairoPDF(filename, width = 10, height = 10)
  print(p)
  dev.off()
  connectivity.zscore
}

#MDS function 
MDSPlot <- function(filename, dataset, targetset, colorscheme = "none", variablename) {
  dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
  target.data <- data.frame(targetset$Slide, factor(targetset[[variablename]]))
  colnames(target.data) <- c("Slide.ID", variablename)
  colnames(dataset.plot) <- c("Slide.ID", "Component.1", "Component.2")
  dataset.plot <- merge(dataset.plot, target.data)
  p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (colorscheme != "none") {
    p <- p + scale_color_manual(values = colorscheme) 
  }
  p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
  CairoPDF(file = filename, height = 7, width = 7)
  print(p)
  dev.off()
}

#Make Excel spreadsheet
DEWorkbook <- function(de.table, filename) { 
  pval.cols <- colnames(de.table) %>% str_detect("P.Value") %>% which
  adj.pval.cols <- colnames(de.table) %>% str_detect("adj.P.Val") %>% which
  coef.cols <- colnames(de.table) %>% str_detect("logFC") %>% which
  
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
  writeDataTable(wb = wb, sheet = 1, x = de.table)
  sig.pvalues <- createStyle(fontColour = "red")
  conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = "<0.005", style = sig.pvalues)
  conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(de.table), rule = "<0.05", style = sig.pvalues)
  conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
  setColWidths(wb, 1, cols = 1, widths = "auto")
  setColWidths(wb, 1, cols = 2, widths = 45)
  setColWidths(wb, 1, cols = 3:12, widths = 15)
  pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
  freezePane(wb, 1, firstRow = TRUE)
  showGridLines(wb, 1, showGridLines = TRUE)
  modifyBaseFont(wb, fontSize = 11)
  saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Plot number of genes at each threshold
DecidePlot <- function(decide.plot, filename, width.plot = 6, height.plot = 7) {
  decide.ggplot <- ggplot() + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")   
  decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = max(Num.Genes) + 110, label = Num.Genes), hjust = -0.3, position = position_dodge(width = 1))
  decide.ggplot <- decide.ggplot + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge") 
  decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = min(Num.Genes) - 110, label = abs(Num.Genes)), hjust = 1.3, position = position_dodge(width = 1))
  decide.ggplot <- decide.ggplot + facet_grid(Test + Num ~ .) 
  decide.ggplot <- decide.ggplot + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  decide.ggplot <- decide.ggplot + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))
  CairoPDF(filename, width = width.plot, height = height.plot)
  print(decide.ggplot)
  dev.off()
}


#Volcano plot ## 
VolcanoPlot <- function(top.table, filename, pval.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue")
{
  top.table$Log.Pvalue <- -log10(top.table[[pval.column]])
  p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = significant)) + xlim(-1.3, 1.3) + ylim(0, 33) + scale_colour_manual(values = c("pink", "lightgrey", "black", "red"))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + xlab(xlabel) + ylab(ylabel)
  CairoPNG(str_c(filename, ".png"), width = 300, height = 300)
  print(p)
  dev.off()
}

GeneBoxplot <- function(lumi.object, gene.symbol) {
  gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
  gene.df <- data.frame(PsychDiagnosis = lumi.object$PsychDiagnosis, Expression = gene.expr, Patient.ID = lumi.object$ID)
  p <- ggplot(gene.df, aes(x = PsychDiagnosis, y = Expression, color = PsychDiagnosis, label = Patient.ID)) + geom_boxplot() + geom_text() + theme_bw()
  p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + ggtitle(gene.symbol) + theme(plot.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
  CairoPDF(gene.symbol, width = 10, height = 10, bg = "transparent")
  print(p)
  dev.off()
}


#GO functions
EnrichrSubmit <- function(dataset, enrichr.terms, colname) {
  dir.create(file.path("./enrichr", colname), showWarnings = FALSE)
  enrichr.data <- map(enrichr.terms, GetEnrichrData, dataset, FALSE)
  enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
  enrichr.data <- enrichr.data[!is.na(enrichr.data)]
  
  names(enrichr.data) <- enrichr.names
  
  map(names(enrichr.data), EnrichrWorkbook, enrichr.data, colname)
  enrichr.data
}

EnrichrWorkbook <- function(database, full.df, colname) {
  dataset <- full.df[[database]]
  
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
  writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
  freezePane(wb, 1, firstRow = TRUE)
  setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
  setColWidths(wb, 1, cols = 2, widths = 45)
  
  dir.create(file.path("./enrichr", colname), recursive = TRUE)
  filename = str_c(file.path("./enrichr", colname, database), ".xlsx")
  saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.height = 5, plot.width = 8) {
  enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
  enrichr.df$Log.pvalue <- -(log10(enrichr.df$Adj.P.value))
  enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
  colnames(enrichr.updown) <- c("Up", "Down")
  enrichr.df <- cbind(enrichr.df, enrichr.updown)
  enrichr.df$Log.Up <- enrichr.df$Log.pvalue * enrichr.df$Up / enrichr.df$Gene.Count
  enrichr.df$Log.Down <- enrichr.df$Log.pvalue * enrichr.df$Down / enrichr.df$Gene.Count
  enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
  enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
  enrichr.df %<>% arrange(Log.pvalue)
  enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
  enrichr.df.plot <- dplyr::select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 
  
  p <- ggplot(select.df.plot, aes(Format.Name, Length, fill = Direction)) + geom_bar(stat = "identity") + geom_text(label = c(as.character(select.df$Format.Name), rep("", nrow(select.df))), hjust = "left", aes(y = 0.1)) + scale_fill_discrete(name = "Direction", labels = c("Up", "Down"))
  p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
  CairoPDF(filename, height = plot.height, width = plot.width)
  print(p)
  dev.off()
}

UpDown <- function(filter.vector, select.df) {
  split.vector <- str_split(filter.vector, ",")[[1]]
  select.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
  enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$z.value) == 1)), "Down" = length(which(sign(enrichr.filter$z.value) == -1)))
  enrichr.vector
}

FilterEnrichr <- function(enrichr.df, size = 200) {
  enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
  enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
  enrichr.filter %<>% dplyr::slice(1:size)
  enrichr.filter
}


###############################################Create lumi object using common NuIDs only ##############################################
#lumibatch1to7.common.NuIDs <- lumiR("./AllRawData/batch1to7AllRaw_common.nuIDs.092319.tsv",lib.mapping = "lumiHumanIDMapping", columnNameGrepPattern = list(exprs='AVG_SIGNAL', se.exprs='BEAD_STDERR', detection='Detection Pval', beadNum='Avg_NBEADS'), checkDupId = TRUE, convertNuID = TRUE, QC = FALSE)
lumibatch1to7.common.NuIDs <- readRDS("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/lumibatch1to7.common.NuIDs.goodrin.492samples.092319.rda")


#Check lumi object info
str(lumibatch1to7.common.NuIDs)

#slideIDs <- read.xlsx("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/transcriptomic_subjects_selection.xlsx", sheet="singletimepoint_noparents")
slideIDs <- read.xlsx("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/transcriptomic_subjects_selection.xlsx", sheet="rm_batch5_singletimepoint_nopar")

selectedsubs_expset.lumi <- lumibatch1to7.common.NuIDs[,lumibatch1to7.common.NuIDs$Slide %in% slideIDs$Slide]
str(selectedsubs_expset.lumi)

missingID <- slideIDs$Slide[!(slideIDs$Slide %in% lumibatch1to7.common.NuIDs$Slide)] 

phenoinfo <- pData(selectedsubs_expset.lumi)

#check age differences
age_group <- lm(Age.at.Collection.Years ~ Diagnosis22q, phenoinfo)
anova(age_group)

phenoinfo %>%
  group_by(Diagnosis22q) %>%
  summarise_at(vars(Age.at.Collection.Years), funs(mean,sd), na.rm=TRUE)

table(phenoinfo$Diagnosis22q)
table(phenoinfo$Diagnosis22q, phenoinfo$Sex2)

#check sex differences
chisq.test(factor(phenoinfo$Sex2), factor(phenoinfo$Diagnosis22q))

phenoinfo %>% dplyr::count(Diagnosis22q, sex) %>% dplyr::group_by(Diagnosis22q) %>% dplyr::mutate(proportion = n/sum(n))

#### Start Pre-processing ####

#START PRE-PROCESSING DATA
lumi.raw <- selectedsubs_expset.lumi
lumi.vst <- lumiT(lumi.raw) #Perform variance stabilized transformation to normalize variance across mean levels of expressed genes
lumi.norm <- lumiN(lumi.vst, method = "rsn") #Normalize with robust spline regression
lumi.cutoff <- detectionCall(lumi.norm) #Get the count of probes which passed the detection threshold per sample
lumi.expr <- lumi.norm[which(lumi.cutoff > 0),] #Subset/Drop any probe where none of the samples passed detection threshold
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated

#Regenerate plots
Histogram("histogram_norm", lumi.expr.annot)
mds.norm <- exprs(lumi.expr.annot) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_diagnosis_norm", mds.norm, pData(lumi.expr.annot), "none", "Diagnosis22q") #label PCs by status

# Remove outliers
connectivity.3group <- ConnectivityPlot("connectivity", lumi.expr.annot, "")
connectivity.outlier.3group <- names(connectivity.3group[connectivity.3group < -3])
remove.key.3group <- !(sampleNames(lumi.expr.annot) %in% connectivity.outlier.3group)
lumi.rmout.3group <- lumi.expr.annot[,remove.key.3group]

final_phenoinfo <- pData(lumi.rmout.3group)

#check age differences
age_group <- lm(Age.at.Collection.Years ~ Diagnosis22q, final_phenoinfo)
anova(age_group)

final_phenoinfo %>%
  group_by(Diagnosis22q) %>%
  summarise_at(vars(Age.at.Collection.Years), funs(mean,sd), na.rm=TRUE)

table(final_phenoinfo$Diagnosis22q)
table(final_phenoinfo$Diagnosis22q, final_phenoinfo$Sex2)

#check sex differences
chisq.test(factor(final_phenoinfo$Sex2), factor(final_phenoinfo$Diagnosis22q))

final_phenoinfo %>% dplyr::count(Diagnosis22q, sex) %>% dplyr::group_by(Diagnosis22q) %>% dplyr::mutate(proportion = n/sum(n))



#Collapse the data by symbol
rmout.collapse.expr.3group <- collapseRows(exprs(lumi.rmout.3group), getSYMBOL(featureNames(lumi.rmout.3group), 'lumiHumanAll.db'), rownames(lumi.rmout.3group)) #collapseRows by symbol
rmout.collapse.3group <- ExpressionSet(assayData = rmout.collapse.expr.3group$datETcollapsed, phenoData = phenoData(lumi.rmout.3group))
lumi.rmout.3group$Sex2 %<>% factor
lumi.rmout.3group$Diagnosis22q %<>% factor

###Remove batch effect AND CELL_TYPE COMPOSITION simultaneously!!
table(phenoinfo$Diagnosis22q, phenoinfo$Batch)
#model.combat.3group <- model.matrix(~ factor(Diagnosis22q) + factor(Sex2) + Age.at.Collection.Years, pData(lumi.expr.annot))

# include batch and cell-type composition using "remove batch effect" function from limma, age and sex can be modeled later on
#model.combat.3group <- model.matrix(~ factor(Diagnosis22q) + factor(Sex2) + Age.at.Collection.Years, pData(lumi.expr.annot))

cibersort_data <- readxl::read_excel("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/regressed_pipeline/rm_batch5_subs/CIBERSORT.Output_Job2.xlsx")
colnames(cibersort_data)[1] <- "Slide" 

pData(rmout.collapse.3group)$Sex2 %<>% as.factor
pheno.for.cibersort <- left_join(pData(rmout.collapse.3group), cibersort_data)
## check "Slide", "Batch", "ID", "Diagnosis22q", "Sex2", "Age.at.Collection.Years") nonmatching

model.design.celltypes <- model.matrix(~ `B cells naive` + `B cells memory` + `Plasma cells` + `T cells CD8` + `T cells CD4 naive` + 
                                         `T cells CD4 memory resting` + `T cells CD4 memory activated` + `T cells follicular helper` +
                                         `T cells regulatory (Tregs)` + `T cells gamma delta` + `NK cells resting` + `NK cells activated` + `Monocytes` +
                                         `Macrophages M0` + `Macrophages M1` + `Macrophages M2` + `Dendritic cells activated` + `Mast cells resting` +
                                         `Mast cells activated`, pheno.for.cibersort) #Make covariate matrix for limma


# Key step to remove batch effect and cell-type abundance at the same time, 
# combat doesn't support continuous variables, equalizes the variance across batches one gene at a time using empirical bayes, better for categorical variables
# because we are correcting for cell-type which is continuous, we can't use combat
rmcelltype.3group <- removeBatchEffect(rmout.collapse.3group, batch = rmout.collapse.3group$Batch, covariates = model.design.celltypes[,-1])
lumi.rmcelltype.3group <- rmout.collapse.3group
exprs(lumi.rmcelltype.3group) <- rmcelltype.3group

# lumi.rmcelltype.3group is new cell-type and batch-regressed expression dataset

pheno.for.cibersort %>%
  group_by(Diagnosis22q) %>%
  summarise_at(vars(Age.at.Collection.Years), funs(mean,sd), na.rm=TRUE)

table(pheno.for.cibersort$Diagnosis22q)
table(pheno.for.cibersort$Diagnosis22q, pheno.for.cibersort$Sex2)

#check sex differences
chisq.test(factor(pheno.for.cibersort$Sex2), factor(pheno.for.cibersort$Diagnosis22q))
pheno.for.cibersort %>% dplyr::count(Diagnosis22q, Sex2) %>% dplyr::group_by(Diagnosis22q) %>% dplyr::mutate(proportion = n/sum(n))

#check batch and RIN
batch.all <- lm(as.integer(Batch) ~ Diagnosis22q, pData(lumi.rmcelltype.3group)) %>% anova %>% tidy
RIN.all <- lm(RIN2 ~ Diagnosis22q, pData(lumi.rmcelltype.3group)) %>% anova %>% tidy
RIN.aov <- aov(as.integer(RIN2) ~ Diagnosis22q, pData(lumi.rmcelltype.3group)) %>% TukeyHSD


pheno.for.cibersort %>%
  group_by(Diagnosis22q) %>%
  summarise_at(vars(RIN2), funs(mean,sd), na.rm=TRUE)

table(pheno.for.cibersort$Diagnosis22q, pheno.for.cibersort$Batch)

#check age differences
age.all <- lm(Age.at.Collection.Years ~ Diagnosis22q, pheno.for.cibersort) %>% anova %>% tidy

pheno.for.cibersort %>%
  group_by(Diagnosis22q) %>%
  summarise_at(vars(Age.at.Collection.Years), funs(mean,sd), na.rm=TRUE)


p <- ggplot(pData(lumi.rmcelltype.3group), aes(x = Diagnosis22q, y = Age.at.Collection.Years, fill = Diagnosis22q)) + geom_violin() + geom_boxplot(width = 0.1) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", signif(age.all$p.value[1], 3), ")", sep = ""))
CairoPDF("./age_boxplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmcelltype.3group), aes(x = Diagnosis22q, y = RIN2)) + geom_boxplot(width = 0.25) + geom_violin()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", round(RIN.all$p.value[1], 3), ")", sep = ""))
CairoPDF("./rin_boxplot.pdf", height = 3, width = 3)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmcelltype.3group), aes(x = Diagnosis22q, fill = factor(Sex2))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle("Sex p = 0.90") 
CairoPDF("./sex_barplot.pdf", height = 3, width = 4)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmcelltype.3group), aes(x = Diagnosis22q, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.position = "none") 
p <- p + ggtitle(paste("Batch (p < ", round(batch.all$p.value[1], 3), ")", sep = "")) 
CairoPDF("./batch_barplot.pdf", height = 3, width = 3)
plot(p)
dev.off()

################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS ###### 
################################################

pheno.for.cibersort$Diagnosis22q <- factor(pheno.for.cibersort$Diagnosis22q, levels = c("control", "22qDup", "22qDel"))

model.design.con.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + as.factor(Sex2) + as.numeric(RIN2)+ Diagnosis22q, pheno.for.cibersort) 
colnames(model.design.con.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "Cont-22qDup", "Cont-22qDel")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes <- lmFit(lumi.rmcelltype.3group, model.design.con.as.ref)
fitb.celltypes <- eBayes(fit.celltypes)

# Controls vs Dup
toptable.ContvsDup.celltypes <- topTable(fitb.celltypes, coef = "Cont-22qDup", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.ContvsDup.celltypes, "./toptable.ContvsDup.celltypes.052020.rda")

# Controls vs Del
toptable.ContvsDel.celltypes <- topTable(fitb.celltypes, coef = "Cont-22qDel", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.ContvsDel.celltypes, "./toptable.ContvsDel.celltypes.052020.rda")

#Make excel datasets

ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description', 'chromosome_name'), filters = 'hgnc_symbol', values = as.character(toptable.ContvsDup.celltypes$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition", "Chromosome")
bm.table %<>% filter(!grepl("CHR", Chromosome))

toptable.ContvsDup.annot.celltypes <- left_join(toptable.ContvsDup.celltypes, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.ContvsDup.annot.celltypes <- toptable.ContvsDup.annot.celltypes[!grepl("PATCH", toptable.ContvsDup.annot.celltypes$Chromosome),]
total_DE_genes <- sum(toptable.ContvsDup.annot.celltypes$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.ContvsDup.annot.celltypes$adj.P.Val <= .05 & toptable.ContvsDup.annot.celltypes$logFC > 0)
total_DE_downgenes <- sum(toptable.ContvsDup.annot.celltypes$adj.P.Val <= .05 & toptable.ContvsDup.annot.celltypes$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.ContvsDup.annot.celltypes, "toptable.ContvsDup_annotated_celltypes_factoredcontrolvars_052020.xlsx")

toptable.ContvsDel.annot.celltypes <- left_join(toptable.ContvsDel.celltypes, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.ContvsDel.annot.celltypes <- toptable.ContvsDel.annot.celltypes[!grepl("PATCH", toptable.ContvsDel.annot.celltypes$Chromosome),]
total_DE_genes <- sum(toptable.ContvsDel.annot.celltypes$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.ContvsDel.annot.celltypes$adj.P.Val <= .05 & toptable.ContvsDel.annot.celltypes$logFC > 0)
total_DE_downgenes <- sum(toptable.ContvsDel.annot.celltypes$adj.P.Val <= .05 & toptable.ContvsDel.annot.celltypes$logFC < 0)

toptable.ContvsDel.annot.celltypes$Chromosome[21] <- 6

#Save workbook and annotation
DEWorkbook(toptable.ContvsDel.annot.celltypes, "toptable.ContvsDel_annotated_celltypes_factoredcontrolvars_052020.xlsx")


# Dup vs Del ### Set deletion as reference, MUST RE-MAKE MODEL DESIGN WITH NEW REFERENCE GROUP
rm(cibersort_data)
cibersort_data <- readxl::read_excel("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/regressed_pipeline/rm_batch5_subs/CIBERSORT.Output_Job2.xlsx")
colnames(cibersort_data)[1] <- "Slide" 
final_phenoinfo$Sex2 <- as.factor(final_phenoinfo$Sex2)

rm(pheno.for.cibersort)
pheno.for.cibersort <- left_join(pData(rmout.collapse.3group), final_phenoinfo)
pheno.for.cibersort$Diagnosis22q <- as.factor(pheno.for.cibersort$Diagnosis22q)


model.design.del.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + as.factor(Sex2) + as.numeric(RIN2)+ Diagnosis22q, pheno.for.cibersort) 
colnames(model.design.del.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "22qDel-22qDup", "22qDel-Cont")


#Limma fit for del vs dup contract with additional covariates
fit.delvsdup.celltypes <- lmFit(lumi.rmcelltype.3group, model.design.del.as.ref)
fitb.delvsdup.celltypes <- eBayes(fit.delvsdup.celltypes)

toptable.DelvsDup.celltypes <- topTable(fitb.delvsdup.celltypes, coef = "22qDel-22qDup", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.DelvsDup.celltypes, "./toptable.DelvsDup.celltypes.520720.rda")

#Save workbook and annotations
toptable.DelvsDup.annot.celltypes <- left_join(toptable.DelvsDup.celltypes, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)

# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.DelvsDup.celltypes <- toptable.DelvsDup.annot.celltypes[!grepl("PATCH", toptable.DelvsDup.annot.celltypes$Chromosome),]
total_DE_genes <- sum(toptable.DelvsDup.annot.celltypes$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.DelvsDup.annot.celltypes$adj.P.Val <= .05 & toptable.DelvsDup.annot.celltypes$logFC > 0)
total_DE_downgenes <- sum(toptable.DelvsDup.annot.celltypes$adj.P.Val <= .05 & toptable.DelvsDup.annot.celltypes$logFC < 0)

toptable.DelvsDup.annot.celltypes$Chromosome[30] <- 6
DEWorkbook(toptable.DelvsDup.annot.celltypes, "toptable.DelvsDup_annotated_celltypes_factoredcontrolvars_052020.xlsx")


## Make Volcano plots

## UPDATED GENE LIST

genes_for_22q <- read_tsv("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/celltype_regression/22q_gene_list.txt")
#install.packages("HGNChelper")
library(HGNChelper)

genes_for_22q.matrix <- genes_for_22q %>% as.matrix()
genelist.genesin22q <- genes_for_22q.matrix[,1]

genelistgenesin22q.hgncfix <- checkGeneSymbols(genelist.genesin22q, unmapped.as.na = FALSE)
genelistgenesin22q.hgnc.genestofix <- genelistgenesin22q.hgncfix[genelistgenesin22q.hgncfix$Approved == FALSE, ] %>% as.data.frame()
genelistgenesin22q.df <- as.data.frame(genes_for_22q$Gene)
colnames(genelistgenesin22q.df) <- "Gene"

genelistgenesin22q.OK <- genelistgenesin22q.df[!(genelistgenesin22q.df$Gene %in% genelistgenesin22q.hgnc.genestofix$x), ] %>% as.data.frame()
colnames(genelistgenesin22q.OK) <- "Gene"
genelistgenesin22q.OK$HGNCUpdatedSymbol <- genelistgenesin22q.OK$Gene

genelistgenesin22q.tofix <- genelistgenesin22q.df[(genelistgenesin22q.df$Gene %in% genelistgenesin22q.hgnc.genestofix$x), ] %>% as.data.frame()
colnames(genelistgenesin22q.tofix) <- "Gene"
colnames(genelistgenesin22q.hgnc.genestofix) <- c("Gene", "Approved", "HGNCUpdatedSymbol")
genelistgenesin22q.hgnc.genestofix.nodupe <- filter(genelistgenesin22q.hgnc.genestofix, !duplicated(Gene))
genelistgenesin22q.tofix.HGNCUpdatedSymbols <- left_join(genelistgenesin22q.tofix, genelistgenesin22q.hgnc.genestofix.nodupe)

genelistgenesin22q.OK.reduce <- dplyr::select(genelistgenesin22q.OK, Gene, HGNCUpdatedSymbol)
genelistgenesin22q.fixed.reduce <- dplyr::select(genelistgenesin22q.tofix.HGNCUpdatedSymbols, Gene, HGNCUpdatedSymbol)
genelistgenesin22q.final <- rbind(genelistgenesin22q.OK.reduce, genelistgenesin22q.fixed.reduce)

## List of updated 22q gene names: genelistgenesin22q.final

genelistgenesin22q.final$cnvgene <- "22q"
final22q.annot <- dplyr::select(genelistgenesin22q.final, HGNCUpdatedSymbol, cnvgene)
colnames(final22q.annot)[1] <- "Symbol"

#Make volcano plots cont vs del, ## MAKE SURE YOU HAVE UPDATED GENE SYMBOLS
toptable.ContvsDel.annot.celltypes <- left_join(toptable.ContvsDel.annot.celltypes, final22q.annot)
significant <- factor(toptable.ContvsDel.annot.celltypes$adj.P.Val < 0.05) %>%
  forcats::fct_recode(sig = "TRUE", notsig = "FALSE") %>%
  as.character
significant[toptable.ContvsDel.annot.celltypes$adj.P.Val < 0.05 & toptable.ContvsDel.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "sig22q_genes"
significant[toptable.ContvsDel.annot.celltypes$adj.P.Val > 0.05 & toptable.ContvsDel.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "nonsig22q_genes"

toptable.ContvsDel.annot.celltypes$significant <- factor(significant,levels = c("nonsig22q_genes", "notsig", "sig", "sig22q_genes")) 

VolcanoPlot(toptable.ContvsDel.annot.celltypes, "volcano_ContvsDel", pval.column = "P.Value", log.column = "logFC")

#Make volcano plots del vs dup, ## MAKE SURE YOU HAVE UPDATED GENE SYMBOLS
toptable.DelvsDup.annot.celltypes <- left_join(toptable.DelvsDup.annot.celltypes, final22q.annot)
significant <- factor(toptable.DelvsDup.annot.celltypes$adj.P.Val < 0.05) %>%
  forcats::fct_recode(sig = "TRUE", notsig = "FALSE") %>%
  as.character
significant[toptable.DelvsDup.annot.celltypes$adj.P.Val < 0.05 & toptable.DelvsDup.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "sig22q_genes"
significant[toptable.DelvsDup.annot.celltypes$adj.P.Val > 0.05 & toptable.DelvsDup.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "nonsig22q_genes"

toptable.DelvsDup.annot.celltypes$significant <- factor(significant,levels = c("nonsig22q_genes", "notsig", "sig", "sig22q_genes")) 

VolcanoPlot(toptable.DelvsDup.annot.celltypes, "volcano_DelvsDup", pval.column = "P.Value", log.column = "logFC")


#Make volcano plots cont vs dup, ## MAKE SURE YOU HAVE UPDATED GENE SYMBOLS
toptable.ContvsDup.annot.celltypes <- left_join(toptable.ContvsDup.annot.celltypes, final22q.annot)
significant <- factor(toptable.ContvsDup.annot.celltypes$adj.P.Val < 0.05) %>%
  forcats::fct_recode(sig = "TRUE", notsig = "FALSE") %>%
  as.character
significant[toptable.ContvsDup.annot.celltypes$adj.P.Val < 0.05 & toptable.ContvsDup.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "sig22q_genes"
significant[toptable.ContvsDup.annot.celltypes$adj.P.Val > 0.05 & toptable.ContvsDup.annot.celltypes$Symbol %in% genelistgenesin22q.final$HGNCUpdatedSymbol] <- "nonsig22q_genes"

toptable.ContvsDup.annot.celltypes$significant <- factor(significant,levels = c("nonsig22q_genes", "notsig", "sig", "sig22q_genes")) 

VolcanoPlot(toptable.ContvsDup.annot.celltypes, "volcano_ContvsDup", pval.column = "P.Value", log.column = "logFC")


#######################################################################################################
###### WGCNA using lumi.rmcelltype.3group which is batch- and celltype-regressed gene expression ######
#######################################################################################################


ModuleWorkbook <- function(module.table, filename) {
  pval.key <- colnames(module.table) %>% str_detect("pvalue") 
  pval.cols <- which(pval.key)
  MM.key <- colnames(module.table) %>% str_detect("MM.") 
  MM.cols <- which(MM.key & !pval.key)
  
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
  writeDataTable(wb = wb, sheet = 1, x = module.table)
  sig.pvalues <- createStyle(fontColour = "red")
  conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(module.table), rule = "<0.05", style = sig.pvalues)
  conditionalFormatting(wb, 1, cols = MM.cols, rows = 1:nrow(module.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
  setColWidths(wb, 1, cols = 1, widths = "auto")
  setColWidths(wb, 1, cols = 2, widths = 45)
  setColWidths(wb, 1, cols = 3:ncol(module.table), widths = "auto")
  pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
  freezePane(wb, 1, firstRow = TRUE)
  showGridLines(wb, 1, showGridLines = TRUE)
  modifyBaseFont(wb, fontSize = 11)
  saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrWorkbook <- function(database, full.df, colname) {
  dataset <- full.df[[database]]
  
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
  writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
  freezePane(wb, 1, firstRow = TRUE)
  setColWidths(wb, 1, cols = 1, widths = 45)
  setColWidths(wb, 1, cols = c(2:ncol(dataset)), widths = "auto")
  
  dir.create(file.path("./enrichr", colname), recursive = TRUE)
  filename = str_c(file.path("./enrichr", colname, database), ".xlsx")
  saveWorkbook(wb, filename, overwrite = TRUE) 
}


GetKscaled <- function(gene.list, module.membership) {
  filter(module.membership, is.element(Symbol, gene.list)) %>% dplyr::select(Symbol, kscaled) %>% arrange(desc(kscaled))
}


powers <- c(1:20) 
# genes are column, samples are row for everytime expression object is called (ie everytime .exprs object is called)!

## START HERE ## lumi.rmcelltype.3group
lumi.rmcelltype.3group.exprs <- exprs(lumi.rmcelltype.3group)
sft <- pickSoftThreshold(t(lumi.rmcelltype.3group.exprs), powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- sft$fitIndices


#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(aes(yintercept = 0.85))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./celltypes_WGCNA_scaleindependence2_power_dupvsdelvscont.pdf", width = 6, height = 6)
print(p)
dev.off()

# #Continue WGCNA with all test samples
softPower <- 12 #was 5 before ### DAN USES AT LEAST 12 to increase the number of modules creates and then have a more stringent merge threshold
adjacency.expr <- adjacency(t(lumi.rmcelltype.3group.exprs), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
TOM <- TOMsimilarity(adjacency.expr, verbose = 5)

dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

#install.packages("flashClust")
library("flashClust")
geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")

#Identify modules using dynamic tree cutting with hybrid clustering
min.module.size = 30 
dynamic.modules <- cutreeDynamic(dendro = geneTree, cutHeight = 0.995, deepSplit = 3, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)

CairoPDF(file = "./celltypes_gene_dendrogram_and_module_colors_signed", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(t(lumi.rmcelltype.3group.exprs), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")

CairoPDF(file = "./celltypes_module_eigengene_clustering_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(lumi.rmcelltype.3group.exprs), dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3, corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

#Use merged eigengenes 
module.colors <- merged.colors
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
ME.genes <- merged.genes

CairoPDF("./celltypes_eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

# Add MEs to pheno data ## dataset now has module eigenegenes! ##
pheno <- pData(rmout.collapse.3group)
all.equal(colnames(lumi.rmcelltype.3group.exprs), pheno$Slide) #to check slide ID order

ME.genes.forPheno <- ME.genes       
ME.genes.forPheno$Slide <- pheno$Slide
pheno.wME <- left_join(pheno, ME.genes.forPheno)

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, dplyr::select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)

gene.info[,3:ncol(gene.info)] <- signif(gene.info[,3:ncol(gene.info)], 3)

#merge module info with biomart info; use manual Entrez annot file rather than biomart
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org");  
bm.table <- getBM(attributes = c('external_gene_name', 'description', 'chromosome_name', 'band'), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition", "Chromosome", "Band")
gene.info.annot <- left_join(gene.info, bm.table) #%>% select(HGNC, Symbol, EntrezGene, Module:kscaled)

gene.module.membership <- data.frame(bicor(t(lumi.rmcelltype.3group.exprs), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(lumi.rmcelltype.3group.exprs))))
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% arrange(Module, desc(kscaled)) %>% left_join(module.membership.pvalue, by = "Symbol")

## CHECK DUPLICATE GENES
dim(module.membership)
module.membership <- distinct(module.membership, Symbol, .keep_all= TRUE)
ModuleWorkbook(module.membership, "celltype_module_membership_22qdupvsdelvsCont_DS1_062520.xlsx")
#check number of genes in each module
table(module.membership$Module)

#Check for Group Differences in Module Eigengene Expression
#table(gene.info$Module) #check module size
EigengeneANOVA <- function(ME.vector, trait.df) {
  trait.final <- mutate(trait.df, ME = ME.vector) #add eigengene as column to phenodata
  #fit linear model with covariates
  trait.anova <- lm(ME ~ Age.at.Collection.Years + as.factor(Sex2) + RIN2 + Diagnosis22q, trait.final) %>% 
    anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  trait.anova$p.value
}

color.values <- unique(module.colors)
anova.diagnosis <- map(ME.genes, EigengeneANOVA, pheno) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "Diagnosis22q")
anova.diagnosis$Module <- colnames(ME.genes)
write.xlsx(anova.diagnosis, "./celltypesWGCNA.ME.DS1.anova.diagnosis.results.xlsx")


#### posthoc pairwise comparisons
sig_modules <- filter(anova.diagnosis, Diagnosis22q <= .056)
sig_module_names <- sig_modules$Module

posthoc_lm <- function(sig_module_names, data_table) {
  module_tmp = data_table[[sig_module_names]]
  
  model_df <- cbind(module_tmp, dplyr::select(data_table, Diagnosis22q))
  print(model_df)
  
  model_aov <- aov(module_tmp ~ as.factor(Diagnosis22q), model_df)  %>% TukeyHSD
  print(model_aov)
  return(model_aov)
}

posthoc.lm.output <- map(sig_module_names, posthoc_lm, pheno.wME) 
names(posthoc.lm.output) <- sig_module_names
capture.output(posthoc.lm.output,file="posthoc.lm.output.doc")


##  make eigeneplots across groups

MEskyblue.module <- dplyr::select(ME.genes, MEskyblue)
MEskyblue.genes.plot <- mutate(MEskyblue.module, Diagnosis22q = pheno$Diagnosis22q) %>%
  gather(Module.Color, Eigengene, -Diagnosis22q) 
MEskyblue.genes.plot$Module.Color %<>% str_replace("ME", "")

p <- ggplot(MEskyblue.genes.plot, aes(x = Diagnosis22q, y = Eigengene, color = Diagnosis22q)) #intialize ggplot object with data
p <- p + geom_boxplot(width = 0.5) #Make boxplot 
p <- p + geom_jitter(width = 0.2) #Show points
p <- p + theme_bw() #Get rid of weird default colors
p <- p + theme(legend.position = "none") #Hide legend
p <- p + theme(panel.grid.major = element_blank()) #Get rid of large gridlines
p <- p + theme(panel.grid.minor = element_blank()) #Get rid of small gridlines
p <- p + theme(axis.title.x = element_blank()) #Get rid of x-axis title
p <- p + theme(axis.title.y = element_blank()) #Get rid of y-axis title
p <- p + theme(panel.border = element_rect(size = 1, color = "turquoise")) #Make the panel border a thicker, turquoise line
p <- p + theme(axis.text = element_text(size = 20)) #Set axis text (tick labels) to larger font size
#p <- p + scale_color_manual(values = c("#000000","#1c9099")) #Override default color scheme
CairoPNG("MEskyblue_eigengene_plot_byDx_DS1.png", height = 450, width = 500, bg = "transparent")
print(p)
dev.off()

######## Test cell-type regressed DE genes against IQ/CT/SA


# rownames(lumi.rmcelltype.3group.exprs) are genes
# colnames(lumi.rmcelltype.3group.exprs) are slide ID
pheno <- pData(lumi.rmcelltype.3group)
all.equal(colnames(lumi.rmcelltype.3group.exprs), pheno$Slide) #to check slide ID order

brain_behavior_data <- read.xlsx("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/revised_brain_beh_data.xlsx")
colnames(brain_behavior_data)[4] <- "visitNum"
str(brain_behavior_data$visitNum)
pheno$visitNum[179] <- 2

brain_behavior_data$visitNum %<>% as.character
all_pheno_brain_beh <- left_join(pheno, brain_behavior_data, by=c("ID","visitNum"))
colnames(all_pheno_brain_beh)

table(all_pheno_brain_beh$Diagnosis22q)

all_pheno_brain_beh$mean_ct <- rowMeans(all_pheno_brain_beh[,c('lh_MeanThickness_thickness','rh_MeanThickness_thickness')])
all_pheno_brain_beh$total_SA <- all_pheno_brain_beh$lh_area + all_pheno_brain_beh$rh_area
all_pheno_brain_beh_df <- tibble::column_to_rownames(all_pheno_brain_beh, "Slide")

pData(lumi.rmcelltype.3group) <- all_pheno_brain_beh_df
validObject(lumi.rmcelltype.3group) #check if sample order is the same in pData and exprs

genes_of_interest <- c("ATRX", "C22orf39", "CDH6", "CHPF2", "CLDN5", "COMT", "CRKL", "DAPK2",	"DGCR2", "DGCR6", "DGCR8", "EDAR", "ESS2", "FBXL5", "FHL3", "FUT7", "GAA", "GNB1L", "HIST1H2BD","HOXC4", "IFT88", "IGFBP4", "KLHL22", "LZTR1", "MED15", "MIR29B2CHG", "PI4KA", "PI4KAP1", "PTPN6",	"RAB11A", "RANBP1", "RTL10", "RTN4R", "SEPTIN5", "SIGLEC10", "SLC22A18", "SLC25A1", "SNAP29",	"TANGO2", "THAP7", "TMEM191A", "TRMT2A", "TXNRD2", "UBAC2", "UFD1", "ZDHHC8", "ZNF319",	"ZNF74", "ADGRE2", "C15orf39", "IPO8", "PLGRKT", "PRSS33", "PTP4A3", "TRPM4", "USP10")

lumi.rmcelltype.genes <- lumi.rmcelltype.3group[genes_of_interest,]
genes_of_interest_expr <- exprs(lumi.rmcelltype.genes) %>% t %>% as_tibble
all_pheno_brain_beh_df$FSIQ %<>% as.numeric 
all_pheno_brain_beh_df$mean_ct %<>% as.numeric 
all_pheno_brain_beh_df$total_SA %<>% as.numeric 

brain_beh_expr_df <- bind_cols(all_pheno_brain_beh_df, genes_of_interest_expr) %>% 
  pivot_longer(ATRX:USP10, names_to = "Gene", values_to = "Expression") %>%
  pivot_longer(one_of(c("mean_ct", "total_SA", "FSIQ")), names_to = "Measure_Name", values_to = "Measure_Value")

FitMeasureLM <- function(all_pheno_rows, all_pheno_group) {
  lm_fit <- lm(Measure_Value ~ Age.at.Collection.Years + Sex2 + Expression * Diagnosis22q, all_pheno_rows)
  lm_tidy <- anova(lm_fit) %>% tidy %>% filter(term != "Residuals")
  lm_tidy$Gene <- all_pheno_group$Gene
  lm_tidy$Measure_Name <- all_pheno_group$Measure_Name
  lm_tidy
}

all_measure_lms <- group_by(brain_beh_expr_df, Gene, Measure_Name) %>% group_map(FitMeasureLM) %>% bind_rows

write.xlsx(all_measure_lms, "./all_measure_lms.xlsx")


######## Test cell-type regressed WGCNA modules genes against IQ/CT/SA

Module_names <- colnames(pheno.wME)[35:62]
Measure_names <- c("total_SA" ,"FSIQ", "mean_ct")
module_plot_data <- left_join(pheno.wME, all_pheno_brain_beh_df)

for(i in 1:28) {
  
  for(j in 1:3) {
    
    gene_plot <- ggplot(data=module_plot_data, aes(x=eval(parse(text = Module_names[i])), y=eval(parse(text = Measure_names[j])), fill = Diagnosis22q)) +
      geom_point() + 
      geom_smooth(method="glm") +
      labs(y = Measure_names[j], x = Module_names[i]) + 
      scale_color_manual(values = c(Deletion = "springgreen3", Control = "lightcoral", Duplication = "cornflowerblue"))
    
    filename <- str_c(Measure_names[j], Module_names[i], ".pdf")
    ggsave(filename,  gene_plot , device = "pdf", width = 5, height = 5)
    print(gene_plot)
  }
}

Module_names <- colnames(pheno.wME)[35:62]
module_exp <- pheno.wME[,35:62]

brain_beh_expr_df <- bind_cols(all_pheno_brain_beh_df, module_exp) %>% 
  pivot_longer(MEmidnightblue:MEgrey, names_to = "Module", values_to = "Expression") %>%
  pivot_longer(one_of(c("mean_ct", "total_SA", "FSIQ")), names_to = "Measure_Name", values_to = "Measure_Value")

FitMeasureLM <- function(all_pheno_rows, all_pheno_group) {
  lm_fit <- lm(Measure_Value ~ Age.at.Collection.Years + Sex2 + Expression * Diagnosis22q, all_pheno_rows)
  lm_tidy <- anova(lm_fit) %>% tidy %>% filter(term != "Residuals")
  lm_tidy$Module <- all_pheno_group$Module
  lm_tidy$Measure_Name <- all_pheno_group$Measure_Name
  lm_tidy
}
module_all_measure_lms <- group_by(brain_beh_expr_df, Module, Measure_Name) %>% group_map(FitMeasureLM) %>% bind_rows

write.xlsx(module_all_measure_lms, "./module_all_measure_lms.xlsx")

######## Test cell-type compositions against IQ/CT/SA
cibersort_results <- read_excel("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/no_fam_effect/final_analysis_pipeline/nonregressed_pipeline/rm_batch5_subs/CIBERSORT.Output_Job2.xlsx")
cibersort_results <- cibersort_results %>% rename(Slide = `Input Sample`)
merged_cibersort <- left_join(final_phenoinfo, cibersort_results)

merged_cibersort <- merged_cibersort %>% rename(`Input Sample` = Slide)
merged_cibersort <- merged_cibersort %>% rename(AGEyears = Age.at.Collection.Years)
merged_cibersort[179,6] <- 2

brain_behavior_data <- read.xlsx("/Users/amylin/Desktop/Users/Amy/Documents/UCLA/Rotations/Bearden_Lab/Transcriptomic_analyses/revised_brain_beh_data.xlsx")
str(brain_behavior_data$timepoint)
colnames(brain_behavior_data)[4] <- "visitNum"
str(brain_behavior_data$visitNum)

all_pheno <- merge(merged_cibersort, brain_behavior_data, by=c("ID","visitNum"))
colnames(all_pheno)

table(merged_cibersort$Diagnosis22q)
table(all_pheno$Diagnosis22q)

all_pheno$mean_ct <- rowMeans(all_pheno[,c('lh_MeanThickness_thickness','rh_MeanThickness_thickness')])
all_pheno$total_SA <- all_pheno$lh_area + all_pheno$rh_area


celltype_names <- c("B cells naive", "B cells memory", "Plasma cells", "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells activated", "Mast cells resting", "Mast cells activated") 


##### Jen and Gil #####
## There will be 6 analyses in total to test the effect of ASD or SCZ on celltype-adjusted DE genes, celltype-adjusted WGCNA modules, and cell-type proportions themselves. 
######################################################################################################################################################################################
## adjusted DE genes (56)           ASD * CNV + ASD + CNV + age + sex + medication1 +  medication2 + etc. <- only including dels and dups (exclude 3's)
## cell-type proportion (19)    ~    
## adjusted WGCNA MEs (28)          SCZ + age + sex +  medication1 +  medication2 + etc. <- only including dels 
######################################################################################################################################################################################

# take first instance of asd or psychosos per ID
# Sort your dataframe by subject ID + ASD status (with ASD status in descending order so 1's occur before 0s), THEN, take first instance per subject in dataframe firsttimept.allsubs <- split(targets.common.NuIDs.pheno, factor(targets.common.NuIDs.pheno$ID)) %>% map_df(slice, 1)
## relevant objects below that were created earlier in the script:
## 
## Adjusted DE genes:
## genes_of_interest <- c("ATRX", "C22orf39", "CDH6", "CHPF2", "CLDN5", "COMT", "CRKL", "DAPK2",	"DGCR2", "DGCR6", "DGCR8", "EDAR", "ESS2", "FBXL5", "FHL3", "FUT7", "GAA", "GNB1L", "HIST1H2BD","HOXC4", "IFT88", "IGFBP4", "KLHL22", "LZTR1", "MED15", "MIR29B2CHG", "PI4KA", "PI4KAP1", "PTPN6",	"RAB11A", "RANBP1", "RTL10", "RTN4R", "SEPTIN5", "SIGLEC10", "SLC22A18", "SLC25A1", "SNAP29",	"TANGO2", "THAP7", "TMEM191A", "TRMT2A", "TXNRD2", "UBAC2", "UFD1", "ZDHHC8", "ZNF319",	"ZNF74", "ADGRE2", "C15orf39", "IPO8", "PLGRKT", "PRSS33", "PTP4A3", "TRPM4", "USP10")
## lumi.rmcelltype.genes <- lumi.rmcelltype.3group[genes_of_interest,] #contains ME data by SlideID with pheno info
## genes_of_interest_expr <- exprs(lumi.rmcelltype.genes) %>% t %>% as_tibble
## 
## Adjusted WGCNA modules:
## Module_names <- colnames(pheno.wME)[35:62]
## pheno.wME #contains ME data by SlideID with pheno info
## 
## Cell-type proportion names:
## celltype_names #contains celltype names to analyse
## merged_cibersort #contains celltype compositions with pheno info


#########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del Psy  No Psy ####### 
#########################################################################

setwd("/Users/ghoftman/Desktop/22qBT")
getwd()


#pData(rmout.collapse.3group)$Sex2 %<>% as.factor
#pheno.for.cibersort <- left_join(pData(rmout.collapse.3group), cibersort_data)
#pheno.22qDel <- pheno[pheno$Diagnosis22q == "22qDel",]

lumi.rmcelltype.del <- lumi.rmcelltype.3group[,lumi.rmcelltype.3group$Diagnosis22q == "22qDel"]
lumi.rmcelltype.del.pheno <- pData(lumi.rmcelltype.del)

med_columns <- read.xlsx("./subject_list_med_cols_forAmy_11092020.xlsx")
med_col_only <- select(med_columns, Slide, SUBJECTID, AD, AED, BZ, AP, ST, NONE)
colnames(med_col_only)[colnames(med_col_only) == 'SUBJECTID'] <- 'ID'
lumi.rmcelltype.del.pheno.med <- left_join(lumi.rmcelltype.del.pheno, med_col_only)
lumi.rmcelltype.del.pheno.med$PSYCHOSIS.DIAGNOS.ANY.TIME <- factor(lumi.rmcelltype.del.pheno.med$PSYCHOSIS.DIAGNOS.ANY.TIME, levels = c("0", "1"))


model.design.NoPsy.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + PSYCHOSIS.DIAGNOS.ANY.TIME, lumi.rmcelltype.del.pheno.med) 
colnames(model.design.NoPsy.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "22qDel.NoPsy-22qDel.Psy")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.NoPsy.Psy <- lmFit(lumi.rmcelltype.del, model.design.NoPsy.as.ref)
fitb.celltypes.NoPsy.Psy <- eBayes(fit.celltypes.NoPsy.Psy)

# 22qDel No Psy vs Psy
toptable.NoPsy.Psy.celltypes <- topTable(fitb.celltypes.NoPsy.Psy, coef = "22qDel.NoPsy-22qDel.Psy", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.NoPsy.Psy.celltypes, "./toptable.NoPsy.Psy.celltypes.121520.rda")

toptable.NoPsy.Psy.annot.celltypes <- left_join(toptable.NoPsy.Psy.celltypes, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.NoPsy.Psy.annot.celltypes <- toptable.NoPsy.Psy.annot.celltypes[!grepl("PATCH", toptable.NoPsy.Psy.annot.celltypes$Chromosome),]
total_DE_genes <- sum(toptable.NoPsy.Psy.annot.celltypes$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.NoPsy.Psy.annot.celltypes$adj.P.Val <= .05 & toptable.NoPsy.Psy.annot.celltypes$logFC > 0)
total_DE_downgenes <- sum(toptable.NoPsy.Psy.annot.celltypes$adj.P.Val <= .05 & toptable.NoPsy.Psy.annot.celltypes$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.NoPsy.Psy.annot.celltypes, "toptable.NoPsy.Psy.annot.celltypes_121520.xlsx")

###### Add Med Columns as Covars into the Model########

#Med columns as Factors
lumi.rmcelltype.del.pheno.med$AD <- factor(lumi.rmcelltype.del.pheno.med$AD, levels = c("0", "1"))
lumi.rmcelltype.del.pheno.med$AED <- factor(lumi.rmcelltype.del.pheno.med$AED, levels = c("0", "1"))
lumi.rmcelltype.del.pheno.med$BZ <- factor(lumi.rmcelltype.del.pheno.med$BZ, levels = c("0", "1"))
lumi.rmcelltype.del.pheno.med$AP <- factor(lumi.rmcelltype.del.pheno.med$AP, levels = c("0", "1"))
lumi.rmcelltype.del.pheno.med$ST <- factor(lumi.rmcelltype.del.pheno.med$ST, levels = c("0", "1"))


model.design.NoPsy.as.ref.meds <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + AD + AED + BZ + AP + ST + PSYCHOSIS.DIAGNOS.ANY.TIME, lumi.rmcelltype.del.pheno.med) 
colnames(model.design.NoPsy.as.ref.meds) <- c("Intercept", "Age", "Sex", "RIN", "AD", "AED", "BZ", "AP", "ST", "22qDel.NoPsy-22qDel.Psy")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.NoPsy.Psy.meds <- lmFit(lumi.rmcelltype.del, model.design.NoPsy.as.ref.meds)
fitb.celltypes.NoPsy.Psy.meds <- eBayes(fit.celltypes.NoPsy.Psy.meds)

# 22qDel No Psy vs Psy w/ Med covars
toptable.NoPsy.Psy.celltypes.meds <- topTable(fitb.celltypes.NoPsy.Psy.meds, coef = "22qDel.NoPsy-22qDel.Psy", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.NoPsy.Psy.celltypes.meds, "./toptable.NoPsy.Psy.celltypes.meds.121520.rda")

toptable.NoPsy.Psy.annot.celltypes.meds <- left_join(toptable.NoPsy.Psy.celltypes.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.NoPsy.Psy.annot.celltypes.meds <- toptable.NoPsy.Psy.annot.celltypes.meds[!grepl("PATCH", toptable.NoPsy.Psy.annot.celltypes.meds$Chromosome),]
total_DE_genes <- sum(toptable.NoPsy.Psy.annot.celltypes.meds$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.NoPsy.Psy.annot.celltypes.meds$adj.P.Val <= .05 & toptable.NoPsy.Psy.annot.celltypes.meds$logFC > 0)
total_DE_downgenes <- sum(toptable.NoPsy.Psy.annot.celltypes.meds$adj.P.Val <= .05 & toptable.NoPsy.Psy.annot.celltypes.meds$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.NoPsy.Psy.annot.celltypes.meds, "toptable.NoPsy.Psy.annot.celltypes.meds_121520.xlsx")


###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del Dup ASD No ASD ###### 
###########################################################################


lumi.rmcelltype.22qdeldup <- lumi.rmcelltype.3group[,lumi.rmcelltype.3group$Diagnosis22q == "22qDel" | lumi.rmcelltype.3group$Diagnosis22q == "22qDup"]
#lumi.rmcelltype.22qdeldup.pheno <- pData(lumi.rmcelltype.22qdeldup)


lumi.rmcelltype.22qdeldup <- lumi.rmcelltype.22qdeldup[,lumi.rmcelltype.22qdeldup$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmcelltype.22qdeldup$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmcelltype.22qdeldup.noNA <- lumi.rmcelltype.22qdeldup[,!is.na(lumi.rmcelltype.22qdeldup$ASD.DIAGNOS.ANY.TIME)]
lumi.rmcelltype.22qdeldup.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcelltype.22qdeldup.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

lumi.rmcelltype.22qdeldup.noNA.pheno <- pData(lumi.rmcelltype.22qdeldup.noNA)

med_columns <- read.xlsx("./subject_list_med_cols_forAmy_11092020.xlsx")
med_col_only <- select(med_columns, Slide, SUBJECTID, AD, AED, BZ, AP, ST, NONE)
colnames(med_col_only)[colnames(med_col_only) == 'SUBJECTID'] <- 'ID'

rm(lumi.rmcelltype.22qdeldup.noNA.pheno.med)
lumi.rmcelltype.22qdeldup.noNA.pheno.med <- left_join(lumi.rmcelltype.22qdeldup.noNA.pheno, med_col_only)
#lumi.rmcelltype.22qdeldup.pheno.med <- filter(lumi.rmcelltype.22qdeldup.pheno.med, ASD.DIAGNOS.ANY.TIME != "3")
#lumi.rmcelltype.22qdeldup.pheno.med$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcelltype.22qdeldup.pheno.med$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

#all.equal(colnames(lumi.rmcelltype.22qdeldup.noNA), lumi.rmcelltype.22qdeldup.pheno.med$ASD.DIAGNOS.ANY.TIME) #to check slide ID order

model.design.22qdeldup.ASDstatus.NoASD.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + Diagnosis22q*ASD.DIAGNOS.ANY.TIME, lumi.rmcelltype.22qdeldup.noNA.pheno.med) 
colnames(model.design.22qdeldup.ASDstatus.NoASD.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "22qDel-22qDup", "NoASD-ASD", "22qDel-Dup:NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.NoASD.ASD <- lmFit(lumi.rmcelltype.22qdeldup.noNA, model.design.22qdeldup.ASDstatus.NoASD.as.ref)
fitb.celltypes.NoASD.ASD <- eBayes(fit.celltypes.NoASD.ASD)

# 22qDel-Dup Main Effect
toptable.22qDelDup.main <- topTable(fitb.celltypes.NoASD.ASD, coef = "22qDel-22qDup", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.22qDelDup.main, "./toptable.22qDelDup.main.121720.rda")

toptable.22qDelDup.main.annot <- left_join(toptable.22qDelDup.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.22qDelDup.main.annot <- toptable.22qDelDup.main.annot[!grepl("PATCH", toptable.22qDelDup.main.annot$Chromosome),]
total_DE_genes_DelDup <- sum(toptable.22qDelDup.main.annot$adj.P.Val <= .05)
total_DE_upgenes_DelDup <- sum(toptable.22qDelDup.main.annot$adj.P.Val <= .05 & toptable.22qDelDup.main.annot$logFC > 0)
total_DE_downgenes_DelDup <- sum(toptable.22qDelDup.main.annot$adj.P.Val <= .05 & toptable.22qDelDup.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.22qDelDup.main.annot, "toptable.22qDelDup.main.annot_121720.xlsx")

# NoASD-ASD Main Effect
toptable.NoASD.ASD.main <- topTable(fitb.celltypes.NoASD.ASD, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.NoASD.ASD.main, "./toptable.NoASD.ASD.main.121720.rda")

toptable.NoASD.ASD.main.annot <- left_join(toptable.NoASD.ASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.NoASD.ASD.main.annot <- toptable.NoASD.ASD.main.annot[!grepl("PATCH", toptable.NoASD.ASD.main.annot$Chromosome),]
total_DE_genes_ASD <- sum(toptable.NoASD.ASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_ASD <- sum(toptable.NoASD.ASD.main.annot$adj.P.Val <= .05 & toptable.NoASD.ASD.main.annot$logFC > 0)
total_DE_downgenes_ASD <- sum(toptable.NoASD.ASD.main.annot$adj.P.Val <= .05 & toptable.NoASD.ASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.NoASD.ASD.main.annot, "toptable.NoASD.ASD.main.annot_121720.xlsx")


# 22qDel-Dup*NoASD-ASD Interaction
toptable.CNV.ASD.inter <- topTable(fitb.celltypes.NoASD.ASD, coef = "22qDel-Dup:NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.CNV.ASD.inter, "./toptable.CNV.ASD.inter.121720.rda")

toptable.CNV.ASD.inter.annot <- left_join(toptable.CNV.ASD.inter, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.CNV.ASD.inter.annot <- toptable.CNV.ASD.inter.annot[!grepl("PATCH", toptable.CNV.ASD.inter.annot$Chromosome),]
total_DE_genes_CNV.ASD <- sum(toptable.CNV.ASD.inter.annot$adj.P.Val <= .05)
total_DE_upgenes_CNV.ASD <- sum(toptable.CNV.ASD.inter.annot$adj.P.Val <= .05 & toptable.CNV.ASD.inter.annot$logFC > 0)
total_DE_downgenes_CNV.ASD <- sum(toptable.CNV.ASD.inter.annot$adj.P.Val <= .05 & toptable.CNV.ASD.inter.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.CNV.ASD.inter.annot, "toptable.CNV.ASD.inter.annot_121720.xlsx")

############# Add Med Columns as Covars into the Model ########3
#Med columns as Factors
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST, levels = c("0", "1"))


model.design.22qdeldup.ASDstatus.NoASD.as.ref.meds <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + AD + AED + BZ + AP + ST + Diagnosis22q*ASD.DIAGNOS.ANY.TIME, lumi.rmcelltype.22qdeldup.noNA.pheno.med) 
colnames(model.design.22qdeldup.ASDstatus.NoASD.as.ref.meds) <- c("Intercept", "Age", "Sex", "RIN", "AD", "AED", "BZ", "AP", "ST", "22qDel-22qDup", "NoASD-ASD", "22qDel-Dup:NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds <- lmFit(lumi.rmcelltype.22qdeldup.noNA, model.design.22qdeldup.ASDstatus.NoASD.as.ref.meds)
fitb.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds <- eBayes(fit.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds)

# 22qDel-Dup CNV main w/ Med covars
toptable.22qDelDup.main.meds <- topTable(fitb.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds, coef = "22qDel-22qDup", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.22qDelDup.main.meds, "./toptable.22qDelDup.main.celltypes.meds.121720.rda")

toptable.22qDelDup.main.annot.celltypes.meds <- left_join(toptable.22qDelDup.main.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.22qDelDup.main.annot.celltypes.meds <- toptable.22qDelDup.main.annot.celltypes.meds[!grepl("PATCH", toptable.22qDelDup.main.annot.celltypes.meds$Chromosome),]
total_DE_genes <- sum(toptable.22qDelDup.main.annot.celltypes.meds$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.22qDelDup.main.annot.celltypes.meds$adj.P.Val <= .05 & toptable.22qDelDup.main.annot.celltypes.meds$logFC > 0)
total_DE_downgenes <- sum(toptable.22qDelDup.main.annot.celltypes.meds$adj.P.Val <= .05 & toptable.22qDelDup.main.annot.celltypes.meds$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.22qDelDup.main.annot.celltypes.meds, "toptable.22qDelDup.main.annot.celltypes.meds_121720.xlsx")

# NoASD-ASD Main Effect
toptable.NoASD.ASD.main.meds <- topTable(fitb.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.NoASD.ASD.main.meds, "./toptable.NoASD.ASD.main.meds.121720.rda")

toptable.NoASD.ASD.main.meds.annot <- left_join(toptable.NoASD.ASD.main.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.NoASD.ASD.main.meds.annot <- toptable.NoASD.ASD.main.meds.annot[!grepl("PATCH", toptable.NoASD.ASD.main.meds.annot$Chromosome),]
total_DE_genes_ASD <- sum(toptable.NoASD.ASD.main.meds$adj.P.Val <= .05)
total_DE_upgenes_ASD <- sum(toptable.NoASD.ASD.main.meds$adj.P.Val <= .05 & toptable.NoASD.ASD.main.meds$logFC > 0)
total_DE_downgenes_ASD <- sum(toptable.NoASD.ASD.main.meds$adj.P.Val <= .05 & toptable.NoASD.ASD.main.meds$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.NoASD.ASD.main.meds.annot, "toptable.NoASD.ASD.main.meds.annot_121720.xlsx")


# 22qDel-Dup*NoASD-ASD Interaction
toptable.CNV.ASD.inter.meds <- topTable(fitb.celltypes.22qdeldup.ASDstatus.NoASD.as.ref.meds, coef = "22qDel-Dup:NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.CNV.ASD.inter.meds, "./toptable.CNV.ASD.inter.meds_121720.rda")

toptable.CNV.ASD.inter.meds.annot <- left_join(toptable.CNV.ASD.inter.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.CNV.ASD.inter.meds.annot <- toptable.CNV.ASD.inter.meds.annot[!grepl("PATCH", toptable.CNV.ASD.inter.meds.annot$Chromosome),]
total_DE_genes_CNV.ASD <- sum(toptable.CNV.ASD.inter.meds.annot$adj.P.Val <= .05)
total_DE_upgenes_CNV.ASD <- sum(toptable.CNV.ASD.inter.meds.annot$adj.P.Val <= .05 & toptable.CNV.ASD.inter.meds.annot$logFC > 0)
total_DE_downgenes_CNV.ASD <- sum(toptable.CNV.ASD.inter.meds.annot$adj.P.Val <= .05 & toptable.CNV.ASD.inter.meds.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.CNV.ASD.inter.meds.annot, "toptable.CNV.ASD.inter.meds.annot_121720.xlsx")



###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del ASD No ASD ###### 
###########################################################################

#Data in Med Columns as Factor
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST, levels = c("0", "1"))

lumi.rmcelltype.22qdel <- lumi.rmcelltype.3group[,lumi.rmcelltype.3group$Diagnosis22q == "22qDel"]

lumi.rmcelltype.22qdel <- lumi.rmcelltype.22qdel[,lumi.rmcelltype.22qdel$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmcelltype.22qdel$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmcelltype.22qdel.noNA <- lumi.rmcelltype.22qdel[,!is.na(lumi.rmcelltype.22qdel$ASD.DIAGNOS.ANY.TIME)]
lumi.rmcelltype.22qdel.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcelltype.22qdel.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

lumi.rmcelltype.22qdel.ASD.noNA.pheno.med <- lumi.rmcelltype.22qdeldup.noNA.pheno.med[lumi.rmcelltype.22qdeldup.noNA.pheno.med$Diagnosis22q == "22qDel",]

model.design.22qdel.ASDstatus.NoASD.as.ref.meds <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + AD + AED + BZ + AP + ST + ASD.DIAGNOS.ANY.TIME, lumi.rmcelltype.22qdel.ASD.noNA.pheno.med) 
colnames(model.design.22qdel.ASDstatus.NoASD.as.ref.meds) <- c("Intercept", "Age", "Sex", "RIN", "AD", "AED", "BZ", "AP", "ST", "NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.22qdel.ASDstatus.NoASD.as.ref.meds <- lmFit(lumi.rmcelltype.22qdel.noNA, model.design.22qdel.ASDstatus.NoASD.as.ref.meds)
fitb.celltypes.22qdel.ASDstatus.NoASD.as.ref.meds <- eBayes(fit.celltypes.22qdel.ASDstatus.NoASD.as.ref.meds)

# 22qDel ASD main w/ Med covars
toptable.22qDel.ASD.main.meds <- topTable(fitb.celltypes.22qdel.ASDstatus.NoASD.as.ref.meds, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.22qDel.ASD.main.meds, "./toptable.22qDel.ASD.main.meds.121720.rda")

toptable.22qDel.ASD.main.meds.annot <- left_join(toptable.22qDel.ASD.main.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.22qDel.ASD.main.meds.annot <- toptable.22qDel.ASD.main.meds.annot[!grepl("PATCH", toptable.22qDel.ASD.main.meds.annot$Chromosome),]
total_DE_genes <- sum(toptable.22qDel.ASD.main.meds.annot$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.22qDel.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.22qDel.ASD.main.meds.annot$logFC > 0)
total_DE_downgenes <- sum(toptable.22qDel.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.22qDel.ASD.main.meds.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.22qDel.ASD.main.meds.annot, "toptable.22qDel.ASD.main.meds.annot_121720.xlsx")


###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Dup ASD No ASD ###### 
###########################################################################

#Make Med Columns data Factors
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AD, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AED, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$BZ, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$AP, levels = c("0", "1"))
lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST <- factor(lumi.rmcelltype.22qdeldup.noNA.pheno.med$ST, levels = c("0", "1"))

lumi.rmcelltype.22qdup <- lumi.rmcelltype.3group[,lumi.rmcelltype.3group$Diagnosis22q == "22qDup"]

lumi.rmcelltype.22qdup <- lumi.rmcelltype.22qdup[,lumi.rmcelltype.22qdup$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmcelltype.22qdup$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmcelltype.22qdup.noNA <- lumi.rmcelltype.22qdup[,!is.na(lumi.rmcelltype.22qdup$ASD.DIAGNOS.ANY.TIME)]
lumi.rmcelltype.22qdup.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcelltype.22qdup.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

lumi.rmcelltype.22qdup.ASD.noNA.pheno.med <- lumi.rmcelltype.22qdeldup.noNA.pheno.med[lumi.rmcelltype.22qdeldup.noNA.pheno.med$Diagnosis22q == "22qDup",]

model.design.22qdup.ASDstatus.NoASD.as.ref.meds <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + AD + AED + AP + ST + ASD.DIAGNOS.ANY.TIME, lumi.rmcelltype.22qdup.ASD.noNA.pheno.med)
#No subjects with BZ meds in this subgroup
colnames(model.design.22qdup.ASDstatus.NoASD.as.ref.meds) <- c("Intercept", "Age", "Sex", "RIN", "AD", "AED", "AP", "ST", "NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.22qdup.ASDstatus.NoASD.as.ref.meds <- lmFit(lumi.rmcelltype.22qdup.noNA, model.design.22qdup.ASDstatus.NoASD.as.ref.meds)
fitb.celltypes.22qdup.ASDstatus.NoASD.as.ref.meds <- eBayes(fit.celltypes.22qdup.ASDstatus.NoASD.as.ref.meds)

# 22qDup ASD main w/ Med covars
toptable.22qDup.ASD.main.meds <- topTable(fitb.celltypes.22qdup.ASDstatus.NoASD.as.ref.meds, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.22qDup.ASD.main.meds, "./toptable.22qDup.ASD.main.meds_121720.rda")

toptable.22qDup.ASD.main.meds.annot <- left_join(toptable.22qDup.ASD.main.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.22qDup.ASD.main.meds.annot <- toptable.22qDup.ASD.main.meds.annot[!grepl("PATCH", toptable.22qDup.ASD.main.meds.annot$Chromosome),]
total_DE_genes <- sum(toptable.22qDup.ASD.main.meds.annot$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.22qDup.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.22qDup.ASD.main.meds.annot$logFC > 0)
total_DE_downgenes <- sum(toptable.22qDup.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.22qDup.ASD.main.meds.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.22qDup.ASD.main.meds.annot, "toptable.22qDup.ASD.main.meds.annot_121720.xlsx")


#########################################################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 2 way interaction Ctl, 22q11 Del, Dup, ASD No ASD ########## 
########  Regressed out cell composition, batch and meds ################### ##########################
#######################################################################################################
getwd()
setwd("/Users/ghoftman/Desktop/22qBT")
lumi.rmcellbatchmed.3group <- readRDS("lumi.rmcelltype.3group.RDA")
lumi.rmcellbatchmed.3group.old <- readRDS("lumi.rmcelltype.3group.RDA")
lumi.rmcellbatchmed.3group.pheno <- pData(lumi.rmcellbatchmed.3group)
lumi.rmcellbatchmed.3group.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0 <- lumi.rmcellbatchmed.3group.pheno$ASD.DIAGNOS.ANY.TIME
lumi.rmcellbatchmed.3group.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0[lumi.rmcellbatchmed.3group.pheno$Diagnosis22q == "control"] <- "0" #only look at control rows within the new column
pData(lumi.rmcellbatchmed.3group.old)$ASD.DIAGNOS.ANY.TIME.CTL.0 = lumi.rmcellbatchmed.3group.pheno$ASD.DIAGNOS.ANY.TIME
pData(lumi.rmcellbatchmed.3group.old)$ASD.DIAGNOS.ANY.TIME.CTL.0[pData(lumi.rmcellbatchmed.3group.old)$Diagnosis22q == "control"] <- "0" #only look at control rows within the new column

lumi.rmcellbatchmed.3group <- lumi.rmcellbatchmed.3group.old


#Remove the "3"s in ASD.DIAGNOS.ANY.TIME (need to keep controls with 3s)
lumi.rmmed.3group <- lumi.rmcellbatchmed.3group[,lumi.rmcellbatchmed.3group$ASD.DIAGNOS.ANY.TIME.CTL.0 == "0" | lumi.rmcellbatchmed.3group$ASD.DIAGNOS.ANY.TIME.CTL.0 == "1"]
#Remove any "NA"s in ASD.DIAGNOS.ANY.TIME (need to keep controls with NAs)
lumi.rmmed.3group.noNA <- lumi.rmmed.3group[,!is.na(lumi.rmmed.3group$ASD.DIAGNOS.ANY.TIME.CTL.0)]
#Make sure ASD.DIAGNOS.ANY.TIME is a factor
lumi.rmmed.3group.noNA$ASD.DIAGNOS.ANY.TIME.CTL.0 <- factor(lumi.rmmed.3group.noNA$ASD.DIAGNOS.ANY.TIME.CTL.0, levels = c("0", "1"))
lumi.rmmed.3group.noNA.pheno <- pData(lumi.rmmed.3group.noNA)
lumi.rmmed.3group.noNA.pheno$Diagnosis22q <- factor(lumi.rmmed.3group.noNA.pheno$Diagnosis22q, levels = c("control", "22qDel", "22qDup"))
lumi.rmmed.3group.noNA$Diagnosis22q <- factor(lumi.rmmed.3group.noNA$Diagnosis22q, levels = c("control", "22qDel", "22qDup"))
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups <- "control"
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups[lumi.rmmed.3group.noNA.pheno$Diagnosis22q == "22qDel" & lumi.rmmed.3group.noNA.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0 == "0"] <- "22qDel.NoASD"
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups[lumi.rmmed.3group.noNA.pheno$Diagnosis22q == "22qDel" & lumi.rmmed.3group.noNA.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0 == "1"] <- "22qDel.ASD"
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups[lumi.rmmed.3group.noNA.pheno$Diagnosis22q == "22qDup" & lumi.rmmed.3group.noNA.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0 == "0"] <- "22qDup.NoASD"
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups[lumi.rmmed.3group.noNA.pheno$Diagnosis22q == "22qDup" & lumi.rmmed.3group.noNA.pheno$ASD.DIAGNOS.ANY.TIME.CTL.0 == "1"] <- "22qDup.ASD"
lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups <- factor(lumi.rmmed.3group.noNA.pheno$ASD.CNV.Groups, levels = c("control", "22qDel.NoASD", "22qDel.ASD", "22qDup.NoASD", "22qDup.ASD"))

rm(model.design.rmmed.3wayInter.ctl.as.ref)

model.design.rmmed.asd.cnv.groups.ctl.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + ASD.CNV.Groups, lumi.rmmed.3group.noNA.pheno) 
colnames(model.design.rmmed.asd.cnv.groups.ctl.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "controlvDel.NoASD", "controlvDel.ASD", "controlvDup.NoASD", "controlvDup.ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.3group.NoASD.ASD <- lmFit(lumi.rmmed.3group.noNA, model.design.rmmed.asd.cnv.groups.ctl.as.ref)
fitb.celltypes.rmmed.3group.NoASD.ASD <- eBayes(fit.celltypes.rmmed.3group.NoASD.ASD)


fit.celltypes.rmmed.3group.NoASD.ASD <- lmFit(lumi.rmmed.3group.noNA, model.design.rmmed.asd.cnv.groups.ctl.as.ref)
contrast.matrix <- makeContrasts(controlvDel.NoASD, controlvDel.ASD, controlvDup.NoASD, controlvDup.ASD, levels = model.design.rmmed.asd.cnv.groups.ctl.as.ref)
fitb.celltypes.rmmed.3group.NoASD.ASD.omni <- contrasts.fit(fit.celltypes.rmmed.3group.NoASD.ASD, contrast.matrix)
fitb.celltypes.rmmed.3group.NoASD.ASD.omni <- eBayes(fitb.celltypes.rmmed.3group.NoASD.ASD.omni)

fitb.omni.pval <- fitb.celltypes.rmmed.3group.NoASD.ASD.omni$F.p.value %>% as.data.frame()
fitb.omni.pval$gene <- rownames(exprs(lumi.rmmed.3group.noNA))
colnames(fitb.omni.pval)[1] <- "p"
fitb.omni.pval$adj.p <- p.adjust(fitb.omni.pval$p, method = "fdr")
##Are the covars included in this analysis?##
write.xlsx(fitb.omni.pval, "omni.adjpval.pairwise_011321.xlsx")

##Make a version of model design w/o covars and run this again##
model.design.rmmed.asd.cnv.groups.ctl.as.ref.nocovars <- model.matrix(~ ASD.CNV.Groups, lumi.rmmed.3group.noNA.pheno) 
colnames(model.design.rmmed.asd.cnv.groups.ctl.as.ref.nocovars) <- c("Intercept", "controlvDel.NoASD", "controlvDel.ASD", "controlvDup.NoASD", "controlvDup.ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.3group.NoASD.ASD.nocovars <- lmFit(lumi.rmmed.3group.noNA, model.design.rmmed.asd.cnv.groups.ctl.as.ref.nocovars)
fitb.celltypes.rmmed.3group.NoASD.ASD.nocovars <- eBayes(fit.celltypes.rmmed.3group.NoASD.ASD.nocovars)


fit.celltypes.rmmed.3group.NoASD.ASD.nocovars <- lmFit(lumi.rmmed.3group.noNA, model.design.rmmed.asd.cnv.groups.ctl.as.ref.nocovars)
contrast.matrix <- makeContrasts(controlvDel.NoASD, controlvDel.ASD, controlvDup.NoASD, controlvDup.ASD, levels = model.design.rmmed.asd.cnv.groups.ctl.as.ref.nocovars)
fitb.celltypes.rmmed.3group.NoASD.ASD.nocovars.omni <- contrasts.fit(fit.celltypes.rmmed.3group.NoASD.ASD.nocovars, contrast.matrix)
fitb.celltypes.rmmed.3group.NoASD.ASD.nocovars.omni <- eBayes(fitb.celltypes.rmmed.3group.NoASD.ASD.nocovars.omni)

fitb.omni.nocovars.pval <- fitb.celltypes.rmmed.3group.NoASD.ASD.nocovars.omni$F.p.value %>% as.data.frame()
fitb.omni.nocovars.pval$gene <- rownames(exprs(lumi.rmmed.3group.noNA))
colnames(fitb.omni.nocovars.pval)[1] <- "p"
fitb.omni.nocovars.pval$adj.p <- p.adjust(fitb.omni.nocovars.pval$p, method = "fdr")
##Are the covars included in this analysis?##
write.xlsx(fitb.omni.nocovars.pval, "omni.adjpval.nocovars.pairwise_011321.xlsx")



#### 5 group Plot for FNBP1####
colnames(lumi.rmmed.3group.noNA.pheno)[11] <- "RNA_Integrity"
lumi.rmmed.3group.noNA.exprs <- exprs(lumi.rmmed.3group.noNA) %>% as.data.frame()
lumi.rmmed.3group.noNA.exprs.t <- t(lumi.rmmed.3group.noNA.exprs) %>% as.data.frame()
lumi.rmmed.3group.noNA.exprs.t$Slide <- rownames(lumi.rmmed.3group.noNA.exprs.t)
lumi.rmmed.3group.noNA.pheno.exprs <- left_join(lumi.rmmed.3group.noNA.pheno, lumi.rmmed.3group.noNA.exprs.t)
lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups <- factor(lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups, levels = c( "22qDel.NoASD", "22qDel.ASD", "control", "22qDup.NoASD", "22qDup.ASD"))

lumi.rmmed.3group.noNA.pheno.exprs.noctl <- lumi.rmmed.3group.noNA.pheno.exprs[lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDel.NoASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDel.ASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDup.NoASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDup.ASD",]

p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs, aes(x = ASD.CNV.Groups, y = FNBP1, color = ASD.CNV.Groups)) + #intialize ggplot object with data
#geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("FNBP1 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=16, angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, angle = 0),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle = 90),
        legend.position = "none") + #Hide legend
        scale_color_manual(values = c("red", "red4", "blue", "green", "green4")) #Override default color scheme
        #legend.title=element_text(size=30),
        #legend.text=element_text(size=30),
        #legend.position = "none")
        #legend.margin = margin(30,0,30,30),
        #legend.key.width = unit(2.5, "line"),
        #legend.key.size = unit(2.8, "lines"))
CairoPNG("FNBP1_plot_byCNVASD_011421.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()

####5 group Plot for WDR1####
p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs, aes(x = ASD.CNV.Groups, y = WDR1, color = ASD.CNV.Groups)) + #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("WDR1 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=16, angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, angle = 0),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle = 90),
        legend.position = "none") + #Hide legend
        scale_color_manual(values = c("red", "red4", "blue", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("WDR1_plot_byCNVASD_011421.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()

#### 5 group Plot for PABPC1 #######
p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs, aes(x = ASD.CNV.Groups, y = CA2, color = ASD.CNV.Groups)) + #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("CA2 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=16, angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, angle = 0),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle = 90),
        legend.position = "none") + #Hide legend
  scale_color_manual(values = c("red", "red4", "blue", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("CA2_plot_byCNVASD_011421.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()

#### Plot for Other genes #######
p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs, aes(x = ASD.CNV.Groups, y = DGCR8, color = ASD.CNV.Groups)) + #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("DGCR8 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=16, angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, angle = 0),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle = 90),
        legend.position = "none") + #Hide legend
  scale_color_manual(values = c("red", "red4", "blue", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("DGCR8_plot_byCNVASD_011421.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()

##### 4 group plots ####
lumi.rmmed.3group.noNA.pheno.exprs.noctl <- lumi.rmmed.3group.noNA.pheno.exprs[lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDel.NoASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDel.ASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDup.NoASD" | lumi.rmmed.3group.noNA.pheno.exprs$ASD.CNV.Groups == "22qDup.ASD",]

p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs.noctl, aes(x = ASD.CNV.Groups, y = FNBP1, color = ASD.CNV.Groups)) + #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("FNBP1 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=20, color = "black", angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=20, color = "black", angle = 0),
        axis.title.x = element_text(size=30, color = "black"),
        axis.title.y = element_text(size=30, color = "black", angle = 90),
        legend.position = "none") + #Hide legend
  scale_color_manual(values = c("red", "red4", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("FNBP1_plot_CNVASD_012621.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()

#### 4 group Plot for WDR1 ####
p <- ggplot(lumi.rmmed.3group.noNA.pheno.exprs.noctl, aes(x = ASD.CNV.Groups, y = WDR1, color = ASD.CNV.Groups)) + #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  geom_boxplot(width = 0.5) + #Make boxplot
  geom_jitter(width = 0.2) + #Show points
  ylab("WDR1 Expression") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2, color = "black"),
        #panel.border = element_blank(),
        axis.text.x = element_text(size=20, color = "black", angle = 0), #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=20, color = "black", angle = 0),
        axis.title.x = element_text(size=30, color = "black"),
        axis.title.y = element_text(size=30, color = "black", angle = 90),
        legend.position = "none") + #Hide legend
  scale_color_manual(values = c("red", "red4", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("WDR1_plot_CNVASD_012621.png", height = 600, width = 700, bg = "transparent")
print(p)
dev.off()


# Ctl-22qDel.ASD Main Effect
toptable.rmmed.Ctl.22qDel.ASD.main <- topTable(fitb.celltypes.rmmed.3group.NoASD.ASD, coef = "control-22qDel.ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.Ctl.22qDel.ASD.main, "./toptable.rmmed.Ctl.22qDel.ASD.main_011221.rda")

toptable.rmmed.Ctl.22qDel.ASD.main.annot <- left_join(toptable.rmmed.Ctl.22qDel.ASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.Ctl.22qDel.ASD.main.annot <- toptable.rmmed.Ctl.22qDel.ASD.main.annot[!grepl("PATCH", toptable.rmmed.Ctl.22qDel.ASD.main.annot$Chromosome),]
total_DE_genes_CtlDelASD <- sum(toptable.rmmed.Ctl.22qDel.ASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_CtlDelASD <- sum(toptable.rmmed.Ctl.22qDel.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDel.ASD.main.annot$logFC > 0)
total_DE_downgenes_CtlDelASD <- sum(toptable.rmmed.Ctl.22qDel.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDel.ASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.Ctl.22qDel.ASD.main.annot, "toptable.rmmed.Ctl.22qDel.ASD.main.annot_011221.xlsx")

# Ctl-22qDup.ASD Main Effect
toptable.rmmed.Ctl.22qDup.ASD.main <- topTable(fitb.celltypes.rmmed.3group.NoASD.ASD, coef = "control-22qDup.ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.Ctl.22qDup.ASD.main, "./toptable.rmmed.Ctl.22qDup.ASD.main_011221.rda")

toptable.rmmed.Ctl.22qDup.ASD.main.annot <- left_join(toptable.rmmed.Ctl.22qDup.ASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.Ctl.22qDup.ASD.main.annot <- toptable.rmmed.Ctl.22qDup.ASD.main.annot[!grepl("PATCH", toptable.rmmed.Ctl.22qDup.ASD.main.annot$Chromosome),]
total_DE_genes_CtlDupASD <- sum(toptable.rmmed.Ctl.22qDup.ASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_CtlDupASD <- sum(toptable.rmmed.Ctl.22qDup.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDup.ASD.main.annot$logFC > 0)
total_DE_downgenes_CtlDupASD <- sum(toptable.rmmed.Ctl.22qDup.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDup.ASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.Ctl.22qDup.ASD.main.annot, "toptable.rmmed.Ctl.22qDup.ASD.main.annot_011221.xlsx")

# Ctl-22qDel.NoASD Main Effect
toptable.rmmed.Ctl.22qDel.NoASD.main <- topTable(fitb.celltypes.rmmed.3group.NoASD.ASD, coef = "control-22qDel.NoASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.Ctl.22qDel.NoASD.main, "./toptable.rmmed.Ctl.22qDel.NoASD.main_011221.rda")

toptable.rmmed.Ctl.22qDel.NoASD.main.annot <- left_join(toptable.rmmed.Ctl.22qDel.NoASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.Ctl.22qDel.NoASD.main.annot <- toptable.rmmed.Ctl.22qDel.NoASD.main.annot[!grepl("PATCH", toptable.rmmed.Ctl.22qDel.NoASD.main.annot$Chromosome),]
total_DE_genes_CtlDelNoASD <- sum(toptable.rmmed.Ctl.22qDel.NoASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_CtlDelNoASD <- sum(toptable.rmmed.Ctl.22qDel.NoASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDel.NoASD.main.annot$logFC > 0)
total_DE_downgenes_CtlDelNoASD <- sum(toptable.rmmed.Ctl.22qDel.NoASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDel.NoASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.Ctl.22qDel.NoASD.main.annot, "toptable.rmmed.Ctl.22qDel.NoASD.main.annot_011221.xlsx")


# Ctl-22qDup.NoASD Main Effect
toptable.rmmed.Ctl.22qDup.NoASD.main <- topTable(fitb.celltypes.rmmed.3group.NoASD.ASD, coef = "control-22qDup.NoASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.Ctl.22qDup.NoASD.main, "./toptable.rmmed.Ctl.22qDup.NoASD.main_011221.rda")

toptable.rmmed.Ctl.22qDup.NoASD.main.annot <- left_join(toptable.rmmed.Ctl.22qDup.NoASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.Ctl.22qDup.NoASD.main.annot <- toptable.rmmed.Ctl.22qDup.NoASD.main.annot[!grepl("PATCH", toptable.rmmed.Ctl.22qDup.NoASD.main.annot$Chromosome),]
total_DE_genes_CtlDupNoASD <- sum(toptable.rmmed.Ctl.22qDup.NoASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_CtlDupNoASD <- sum(toptable.rmmed.Ctl.22qDup.NoASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDup.NoASD.main.annot$logFC > 0)
total_DE_downgenes_CtlDupNoASD <- sum(toptable.rmmed.Ctl.22qDup.NoASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.Ctl.22qDup.NoASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.Ctl.22qDup.NoASD.main.annot, "toptable.rmmed.Ctl.22qDup.NoASD.main.annot_011221.xlsx")



#########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del Psy  No Psy ###########
## cell comp, batch and meds regressed out      ####### 
#########################################################################

setwd("/Users/ghoftman/Desktop/22qBT")
getwd()

lumi.rmcellbatchmed.3group <- readRDS("/Users/ghoftman/Desktop/22qBT/lumi.rmcelltype.3group.RDA")

lumi.rmcellbatchmed.del <- lumi.rmcellbatchmed.3group[,lumi.rmcellbatchmed.3group$Diagnosis22q == "22qDel"]
lumi.rmcellbatchmed.del.pheno <- pData(lumi.rmcellbatchmed.del)

table(lumi.rmcellbatchmed.del.pheno$Diagnosis22q, lumi.rmcellbatchmed.del.pheno$PSYCHOSIS.DIAGNOS.ANY.TIME)


med_columns <- read.xlsx("./subject_list_med_cols_forAmy_11092020.xlsx")
med_col_only <- select(med_columns, Slide, SUBJECTID, AD, AED, BZ, AP, ST, NONE)
colnames(med_col_only)[colnames(med_col_only) == 'SUBJECTID'] <- 'ID'
#lumi.rmcelltype.del.pheno.med <- left_join(lumi.rmcelltype.del.pheno, med_col_only)
lumi.rmcellbatchmed.del.pheno$PSYCHOSIS.DIAGNOS.ANY.TIME <- factor(lumi.rmcellbatchmed.del.pheno$PSYCHOSIS.DIAGNOS.ANY.TIME, levels = c("0", "1"))

model.design.rmmed.NoPsy.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + PSYCHOSIS.DIAGNOS.ANY.TIME, lumi.rmcellbatchmed.del.pheno) 
colnames(model.design.rmmed.NoPsy.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "22qDel.NoPsy-22qDel.Psy")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.NoPsy.Psy <- lmFit(lumi.rmcellbatchmed.del, model.design.rmmed.NoPsy.as.ref)
fitb.celltypes.rmmed.NoPsy.Psy <- eBayes(fit.celltypes.rmmed.NoPsy.Psy)

# 22qDel No Psy vs Psy
toptable.rmmed.NoPsy.Psy.celltypes <- topTable(fitb.celltypes.rmmed.NoPsy.Psy, coef = "22qDel.NoPsy-22qDel.Psy", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.NoPsy.Psy.celltypes, "./toptable.rmmed.NoPsy.Psy.celltypes_122220.rda")

toptable.rmmed.NoPsy.Psy.celltypes.annot <- left_join(toptable.rmmed.NoPsy.Psy.celltypes, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.NoPsy.Psy.celltypes.annot <- toptable.rmmed.NoPsy.Psy.celltypes.annot[!grepl("PATCH", toptable.rmmed.NoPsy.Psy.celltypes.annot$Chromosome),]
total_DE_genes <- sum(toptable.rmmed.NoPsy.Psy.celltypes.annot$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.rmmed.NoPsy.Psy.celltypes.annot$adj.P.Val <= .05 & toptable.rmmed.NoPsy.Psy.celltypes.annot$logFC > 0)
total_DE_downgenes <- sum(toptable.rmmed.NoPsy.Psy.celltypes.annot$adj.P.Val <= .05 & toptable.rmmed.NoPsy.Psy.celltypes.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.NoPsy.Psy.celltypes.annot, "toptable.rmmed.NoPsy.Psy.celltypes.annot_122220.xlsx")


###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del Dup ASD No ASD 
########  Regressed out cell composition, batch and meds ################### 
###########################################################################

lumi.rmmed.22qdeldup <- lumi.rmcellbatchmed.3group[,lumi.rmcellbatchmed.3group$Diagnosis22q == "22qDel" | lumi.rmcellbatchmed.3group$Diagnosis22q == "22qDup"]
#lumi.rmcelltype.22qdeldup.pheno <- pData(lumi.rmcelltype.22qdeldup)

lumi.rmmed.22qdeldup <- lumi.rmmed.22qdeldup[,lumi.rmmed.22qdeldup$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmmed.22qdeldup$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmmed.22qdeldup.noNA <- lumi.rmmed.22qdeldup[,!is.na(lumi.rmmed.22qdeldup$ASD.DIAGNOS.ANY.TIME)]
lumi.rmmed.22qdeldup.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmmed.22qdeldup.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))
lumi.rmmed.22qdeldup.noNA.pheno <- pData(lumi.rmmed.22qdeldup.noNA)

model.design.22qdeldup.ASDstatus.rmmed.NoASD.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + Diagnosis22q*ASD.DIAGNOS.ANY.TIME, lumi.rmmed.22qdeldup.noNA.pheno) 
colnames(model.design.22qdeldup.ASDstatus.rmmed.NoASD.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "22qDel-22qDup", "NoASD-ASD", "22qDel-Dup:NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.NoASD.ASD <- lmFit(lumi.rmmed.22qdeldup.noNA, model.design.22qdeldup.ASDstatus.rmmed.NoASD.as.ref)
fitb.celltypes.rmmed.NoASD.ASD <- eBayes(fit.celltypes.rmmed.NoASD.ASD)

# 22qDel-Dup Main Effect
toptable.rmmed.22qDelDup.main <- topTable(fitb.celltypes.rmmed.NoASD.ASD, coef = "22qDel-22qDup", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.22qDelDup.main, "./toptable.rmmed.22qDelDup.main_122220.rda")

toptable.rmmed.22qDelDup.main.annot <- left_join(toptable.rmmed.22qDelDup.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.22qDelDup.main.annot <- toptable.rmmed.22qDelDup.main.annot[!grepl("PATCH", toptable.rmmed.22qDelDup.main.annot$Chromosome),]
total_DE_genes_DelDup <- sum(toptable.rmmed.22qDelDup.main.annot$adj.P.Val <= .05)
total_DE_upgenes_DelDup <- sum(toptable.rmmed.22qDelDup.main.annot$adj.P.Val <= .05 & toptable.rmmed.22qDelDup.main.annot$logFC > 0)
total_DE_downgenes_DelDup <- sum(toptable.rmmed.22qDelDup.main.annot$adj.P.Val <= .05 & toptable.rmmed.22qDelDup.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.22qDelDup.main.annot, "toptable.rmmed.22qDelDup.main.annot_122220.xlsx")

# NoASD-ASD Main Effect
toptable.rmmed.NoASD.ASD.main <- topTable(fitb.celltypes.rmmed.NoASD.ASD, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.NoASD.ASD.main, "./toptable.rmmed.NoASD.ASD.main_122220.rda")

toptable.rmmed.NoASD.ASD.main.annot <- left_join(toptable.rmmed.NoASD.ASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.NoASD.ASD.main.annot <- toptable.rmmed.NoASD.ASD.main.annot[!grepl("PATCH", toptable.rmmed.NoASD.ASD.main.annot$Chromosome),]
total_DE_genes_ASD <- sum(toptable.rmmed.NoASD.ASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes_ASD <- sum(toptable.rmmed.NoASD.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.NoASD.ASD.main.annot$logFC > 0)
total_DE_downgenes_ASD <- sum(toptable.rmmed.NoASD.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.NoASD.ASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.NoASD.ASD.main.annot, "toptable.rmmed.NoASD.ASD.main.annot_122220.xlsx")


# 22qDel-Dup*NoASD-ASD Interaction
toptable.rmmed.CNVxASD <- topTable(fitb.celltypes.rmmed.NoASD.ASD, coef = "22qDel-Dup:NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.CNVxASD, "./toptable.rmmed.CNVxASD_122220.rda")

toptable.rmmed.CNVxASD.annot <- left_join(toptable.rmmed.CNVxASD, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.CNVxASD.annot <- toptable.rmmed.CNVxASD.annot[!grepl("PATCH", toptable.rmmed.CNVxASD.annot$Chromosome),]
total_DE_genes_CNV.ASD <- sum(toptable.rmmed.CNVxASD.annot$adj.P.Val <= .05)
total_DE_upgenes_CNV.ASD <- sum(toptable.rmmed.CNVxASD.annot$adj.P.Val <= .05 & toptable.rmmed.CNVxASD.annot$logFC > 0)
total_DE_downgenes_CNV.ASD <- sum(toptable.rmmed.CNVxASD.annot$adj.P.Val <= .05 & toptable.rmmed.CNVxASD.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.CNVxASD.annot, "toptable.rmmed.CNVxASD.annot_122220.xlsx")


###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Del ASD No ASD ######
######### cell composition, batch and meds regressed out ##############
###########################################################################

lumi.rmcellbatchmed.del <- lumi.rmcellbatchmed.3group[,lumi.rmcellbatchmed.3group$Diagnosis22q == "22qDel"]
lumi.rmcellbatchmed.del.ASD <- lumi.rmcellbatchmed.del[,lumi.rmcellbatchmed.del$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmcellbatchmed.del$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmcellbatchmed.del.ASD.noNA <- lumi.rmcellbatchmed.del.ASD[,!is.na(lumi.rmcellbatchmed.del.ASD$ASD.DIAGNOS.ANY.TIME)]
lumi.rmcellbatchmed.del.ASD.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcellbatchmed.del.ASD.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))
lumi.rmcellbatchmed.del.ASD.noNA.pheno <- pData(lumi.rmcellbatchmed.del.ASD.noNA)

model.design.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + ASD.DIAGNOS.ANY.TIME, lumi.rmcellbatchmed.del.ASD.noNA.pheno) 
colnames(model.design.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds) <- c("Intercept", "Age", "Sex", "RIN", "NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds <- lmFit(lumi.rmcellbatchmed.del.ASD.noNA, model.design.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds)
fitb.celltypes.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds <- eBayes(fit.celltypes.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds)

# 22qDel ASD main w/ Med covars
toptable.rmmed.22qDel.ASD.main.meds <- topTable(fitb.celltypes.rmmed.22qdel.ASDstatus.NoASD.as.ref.meds, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.22qDel.ASD.main.meds, "./toptable.rmmed.22qDel.ASD.main.meds_122220.rda")

toptable.rmmed.22qDel.ASD.main.meds.annot <- left_join(toptable.rmmed.22qDel.ASD.main.meds, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.22qDel.ASD.main.meds.annot <- toptable.rmmed.22qDel.ASD.main.meds.annot[!grepl("PATCH", toptable.rmmed.22qDel.ASD.main.meds.annot$Chromosome),]
total_DE_genes <- sum(toptable.rmmed.22qDel.ASD.main.meds.annot$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.rmmed.22qDel.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.rmmed.22qDel.ASD.main.meds.annot$logFC > 0)
total_DE_downgenes <- sum(toptable.rmmed.22qDel.ASD.main.meds.annot$adj.P.Val <= .05 & toptable.rmmed.22qDel.ASD.main.meds.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.22qDel.ASD.main.meds.annot, "toptable.rmmed.22qDel.ASD.main.meds.annot_122220.xlsx")


###########################################################################
###### DIFFERENTIAL EXPRESSION ANALYSIS for 22q11 Dup ASD No ASD ###### 
#### regressed out cell composition, batch and meds ####################
###########################################################################

lumi.rmcellbatchmed.ASD.dup <- lumi.rmcellbatchmed.3group[,lumi.rmcellbatchmed.3group$Diagnosis22q == "22qDup"]

lumi.rmcellbatchmed.ASD.dup <- lumi.rmcellbatchmed.ASD.dup[,lumi.rmcellbatchmed.ASD.dup$ASD.DIAGNOS.ANY.TIME == "0" | lumi.rmcellbatchmed.ASD.dup$ASD.DIAGNOS.ANY.TIME == "1"]
lumi.rmcellbatchmed.ASD.dup.noNA <- lumi.rmcellbatchmed.ASD.dup[,!is.na(lumi.rmcellbatchmed.ASD.dup$ASD.DIAGNOS.ANY.TIME)]
lumi.rmcellbatchmed.ASD.dup.noNA$ASD.DIAGNOS.ANY.TIME <- factor(lumi.rmcellbatchmed.ASD.dup.noNA$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

lumi.rmcellbatchmed.ASD.dup.noNA.pheno <- pData(lumi.rmcellbatchmed.ASD.dup.noNA)

model.design.rmmed.22qdup.ASDstatus.NoASD.as.ref <- model.matrix(~ as.numeric(Age.at.Collection.Years) + Sex2 + as.numeric(RIN2) + ASD.DIAGNOS.ANY.TIME, lumi.rmcellbatchmed.ASD.dup.noNA.pheno)
colnames(model.design.rmmed.22qdup.ASDstatus.NoASD.as.ref) <- c("Intercept", "Age", "Sex", "RIN", "NoASD-ASD")

#Limma fit on cell-type and batch-corrected gene expression for additional covariates
fit.celltypes.rmmed.22qdup.ASDstatus.NoASD.as.ref <- lmFit(lumi.rmcellbatchmed.ASD.dup.noNA, model.design.rmmed.22qdup.ASDstatus.NoASD.as.ref)
fitb.celltypes.rmmed.22qdup.ASDstatus.NoASD.as.ref <- eBayes(fit.celltypes.rmmed.22qdup.ASDstatus.NoASD.as.ref)

# 22qDup ASD main w/ Med covars
toptable.rmmed.22qDup.ASD.main <- topTable(fitb.celltypes.rmmed.22qdup.ASDstatus.NoASD.as.ref, coef = "NoASD-ASD", n = Inf, sort.by = "p") %>%
  mutate(Symbol = rownames(.)) %>% as_tibble
saveRDS(toptable.rmmed.22qDup.ASD.main, "./toptable.rmmed.22qDup.ASD.main_122220.rda")

toptable.rmmed.22qDup.ASD.main.annot <- left_join(toptable.rmmed.22qDup.ASD.main, bm.table) %>% 
  dplyr::select(Symbol, Definition, Chromosome, logFC, P.Value, adj.P.Val, t, B, AveExpr) %>%
  arrange(P.Value)
# remove genes where chromosome name includes "XXX_PATCH" as they are duplicates
toptable.rmmed.22qDup.ASD.main.annot <- toptable.rmmed.22qDup.ASD.main.annot[!grepl("PATCH", toptable.rmmed.22qDup.ASD.main.annot$Chromosome),]
total_DE_genes <- sum(toptable.rmmed.22qDup.ASD.main.annot$adj.P.Val <= .05)
total_DE_upgenes <- sum(toptable.rmmed.22qDup.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.22qDup.ASD.main.annot$logFC > 0)
total_DE_downgenes <- sum(toptable.rmmed.22qDup.ASD.main.annot$adj.P.Val <= .05 & toptable.rmmed.22qDup.ASD.main.annot$logFC < 0)

#Save workbook and annotations
DEWorkbook(toptable.rmmed.22qDup.ASD.main.annot, "toptable.rmmed.22qDup.ASD.main.annot_122220.xlsx")



######################### WGCNA for 22q11 Del Psy No Psy #########################################
########################   with cell composition, meds and batch regressed out ################

EigengeneANOVA.rmmeds.psy.nopsy <- function(ME.vector, trait.df) {
  trait.final <- mutate(trait.df, ME = ME.vector) #add eigengene as column to phenodata
  #fit linear model with covariates
  trait.anova <- lm(ME ~ Age.at.Collection.Years + as.factor(Sex2) + RIN2 + PSYCHOSIS.DIAGNOS.ANY.TIME, trait.final) %>% 
    anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  trait.anova$p.value
}

EigengeneANOVA.rmmeds.psy.nopsy.nogrey <- function(ME.vector, trait.df) {
  trait.final <- mutate(trait.df, ME = ME.vector) #add eigengene as column to phenodata
  #fit linear model with covariates
  trait.anova <- lm(ME ~ Age.at.Collection.Years + as.factor(Sex2) + RIN2 + PSYCHOSIS.DIAGNOS.ANY.TIME, trait.final) %>% 
    anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  trait.anova$p.value
}

# Add MEs to pheno data ## dataset now has module eigenegenes! #

#pheno <- pData(rmout.collapse.3group)
#all.equal(rownames(ME.genes), pheno$Slide) #to check slide ID order
#ME.genes.forPheno <- ME.genes       
#ME.genes.forPheno$Slide <- pheno$Slide
#pheno.wME <- left_join(pheno, ME.genes.forPheno)

pheno.wME.rmmed.update <- readRDS("pheno.wME.2021.rda")
pheno.wME.rmmed.22qDel.update <- pheno.wME.rmmed.update[pheno.wME.rmmed.update$Diagnosis22q == "22qDel",] 
#pheno.Del <- pheno[pheno$Diagnosis22q == "22qDel",]
ME.genes.Del <- ME.genes[rownames(ME.genes) %in% pheno.Del$Slide,]
all.equal(rownames(ME.genes.Del), pheno.wME.rmmed.22qDel.update$Slide)
pheno.wME.rmmed.22qDel.update$PSYCHOSIS.DIAGNOS.ANY.TIME <- factor(pheno.wME.rmmed.22qDel.update$PSYCHOSIS.DIAGNOS.ANY.TIME, levels = c("0", "1"))

ME.genes.Del.nogrey <- ME.genes.Del[,-28]

### ANOVA for FDR corrected p vals with and w/o the grey catch all module ###
color.values <- unique(module.colors)
anova.diagnosis.rmmed.psy.nopsy <- map(ME.genes.Del, EigengeneANOVA.rmmeds.psy.nopsy, pheno.wME.rmmed.22qDel.update) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.rmmed.psy.nopsy) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel.NoPsy-22qDel.Psy")
anova.diagnosis.rmmed.psy.nopsy$Module <- colnames(ME.genes)
write.xlsx(anova.diagnosis.rmmed.psy.nopsy, "./celltypesWGCNA.ME.anova.22qDel.rmmed.psy.nopsy.update_012621.xlsx")

anova.diagnosis.rmmed.psy.nopsy.nogrey <- map(ME.genes.Del.nogrey, EigengeneANOVA.rmmeds.psy.nopsy.nogrey, pheno.wME.rmmed.22qDel.update) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.rmmed.psy.nopsy.nogrey) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel.NoPsy-22qDel.Psy")
anova.diagnosis.rmmed.psy.nopsy.nogrey$Module <- colnames(ME.genes.Del.nogrey)
write.xlsx(anova.diagnosis.rmmed.psy.nopsy.nogrey, "./celltypesWGCNA.ME.anova.22qDel.rmed.psy.nopsy.nogrey.update_012621.xlsx")

rm(pheno.wME.rmmed.22qDel.update.long)
### Convert to long form for plotting ###
pheno.wME.rmmed.22qDel.update.long <- gather(pheno.wME.rmmed.22qDel.update, Module, Expression, MEyellow:MEcyan, factor_key = TRUE)
pheno.wME.rmmed.22qDel.update.long$Module <- str_replace_all(pheno.wME.rmmed.22qDel.update.long$Module, "ME", "")
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce <-  pheno.wME.rmmed.22qDelDup.ASD.long[!(pheno.wME.rmmed.22qDelDup.ASD.long$Module == "grey"), ]
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce <- left_join(pheno.wME.rmmed.22qDelDup.ASD.long.reduce, Kangmodule.key)
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce$ModuleNumberName <- factor(pheno.wME.rmmed.22qDelDup.ASD.long.reduce$ModuleNumberName, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18"))

#####  Create the PsyGroups column and then set values to names ####
pheno.wME.rmmed.22qDel.update.long$PsyGroups <- "22qDel.NoPsy"
pheno.wME.rmmed.22qDel.update.long$PsyGroups[pheno.wME.rmmed.22qDel.update.long$Diagnosis22q == "22qDel" & pheno.wME.rmmed.22qDel.update.long$PSYCHOSIS.DIAGNOS.ANY.TIME == "1"] <- "22qDel.Psy"
pheno.wME.rmmed.22qDel.update.long$PsyGroups <- factor(pheno.wME.rmmed.22qDel.update.long$PsyGroups, levels = c("22qDel.NoPsy", "22qDel.Psy"))

#### Plot the ME box plots as a facet ###
p <- ggplot(pheno.wME.rmmed.22qDel.update.long, aes(x = PsyGroups, y = Expression, color = PsyGroups)) +
  #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  facet_wrap(~ Module, ncol = 5) +
  geom_boxplot(width = 0.7) + #Make boxplot
  geom_jitter(width = 0.1) + #Show points
  ylab("Module Eigengene") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5, color = "black"),
        panel.spacing = unit(2, "lines"),
        #panel.border = element_blank(),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(size = 16, face = "bold", color = "black", angle = 45, hjust=0.6, vjust=0.6), 
        #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, color = "black", angle = 0),
        axis.title.x = element_text(size=30, color = "black", angle = 0),
        axis.title.y = element_text(size=30, color = "black", angle = 90),
        legend.position = "none") +
  scale_color_manual(values = c("red", "red4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("ME_plot_Psy_facet_012621.png", height = 1200, width = 1200, bg = "transparent")
print(p)
dev.off()


#################### WGCNA for 22q11 Del Dup ASD No ASD #######################################
########################   with cell composition, meds and batch regressed out ################

EigengeneANOVA.deldup.asd.noasd <- function(ME.vector, trait.df) {
  trait.final <- mutate(trait.df, ME = ME.vector) #add eigengene as column to phenodata
  #fit linear model with covariates
  trait.anova <- lm(ME ~ Age.at.Collection.Years + as.factor(Sex2) + RIN2 + Diagnosis22q*ASD.DIAGNOS.ANY.TIME, trait.final) %>% 
    anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  trait.anova$p.value
}

EigengeneANOVA.deldup.asd.noasd.nogrey <- function(ME.vector, trait.df) {
  trait.final <- mutate(trait.df, ME = ME.vector) #add eigengene as column to phenodata
  #fit linear model with covariates
  trait.anova <- lm(ME ~ Age.at.Collection.Years + as.factor(Sex2) + RIN2 + Diagnosis22q*ASD.DIAGNOS.ANY.TIME, trait.final) %>% 
    anova %>% 
    tidy %>%
    filter(!is.na(p.value))
  trait.anova$p.value
}

# Add MEs to pheno data ## dataset now has module eigenegenes! #

#pheno <- pData(rmout.collapse.3group)
#all.equal(rownames(ME.genes), pheno$Slide) #to check slide ID order
#ME.genes.forPheno <- ME.genes       
#ME.genes.forPheno$Slide <- pheno$Slide
#pheno.wME <- left_join(pheno, ME.genes.forPheno)

pheno.wME.rmmed.22qDelDup.update <- pheno.wME.rmmed.update[pheno.wME.rmmed.update$Diagnosis22q == "22qDel" | pheno.wME.rmmed.update$Diagnosis22q == "22qDup",] 
pheno.DelDup <- pheno[pheno$Diagnosis22q == "22qDel" | pheno$Diagnosis22q == "22qDup",]

pheno.wME.rmmed.22qDelDup.update.ASD <- pheno.wME.rmmed.22qDelDup.update[pheno.wME.rmmed.22qDelDup.update$ASD.DIAGNOS.ANY.TIME == "0" | pheno.wME.rmmed.22qDelDup.update$ASD.DIAGNOS.ANY.TIME == "1",]
#pheno.wME.rmmed.22qDelDup.ASD.noNA <- pheno.wME.rmmed.22qDelDup.ASD[!is.na(pheno.wME.rmmed.22qDelDup.ASD$ASD.DIAGNOS.ANY.TIME),]
pheno.wME.rmmed.22qDelDup.update.ASD$ASD.DIAGNOS.ANY.TIME <- factor(pheno.wME.rmmed.22qDelDup.update.ASD$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))


pheno.DelDup.ASD <- pheno.DelDup[pheno.DelDup$ASD.DIAGNOS.ANY.TIME == "0" | pheno.DelDup$ASD.DIAGNOS.ANY.TIME == "1",]
#pheno.DelDup.ASD.noNA <- pheno.DelDup.ASD[,!is.na(pheno.DelDup.ASD$ASD.DIAGNOS.ANY.TIME)]
pheno.DelDup.ASD$ASD.DIAGNOS.ANY.TIME <- factor(pheno.DelDup.ASD$ASD.DIAGNOS.ANY.TIME, levels = c("0", "1"))

ME.genes.DelDup.ASD <- ME.genes[rownames(ME.genes) %in% pheno.DelDup.ASD$Slide,]

ME.genes.DelDup.ASD.nogrey <- ME.genes.DelDup.ASD[,-28]  

all.equal(rownames(ME.genes.DelDup.ASD), pheno.wME.rmmed.22qDelDup.update.ASD$Slide)
all.equal(rownames(ME.genes.DelDup.ASD), pheno.DelDup.ASD$Slide)

### ANOVA FDR corrected p vals with and w/o the grey catch all module ###
color.values <- unique(module.colors)
anova.diagnosis.deldup.asd.noasd.nogrey <- map(ME.genes.DelDup.ASD.nogrey, EigengeneANOVA.deldup.asd.noasd.nogrey, pheno.wME.rmmed.22qDelDup.update.ASD) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.deldup.asd.noasd.nogrey) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel-22qDup", "NoASD-ASD", "CNV:ASD")
anova.diagnosis.deldup.asd.noasd.nogrey$Module <- colnames(ME.genes.DelDup.ASD.nogrey)
write.xlsx(anova.diagnosis.deldup.asd.noasd.nogrey, "./celltypesWGCNA.ME.anova.diagnosis.deldup.asd.noasd.nogrey_012621.xlsx")

anova.diagnosis.deldup.asd.noasd <- map(ME.genes.DelDup.ASD, EigengeneANOVA.deldup.asd.noasd, pheno.wME.rmmed.22qDelDup.update.ASD) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(p.adjust, method = "fdr") %>% 
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.deldup.asd.noasd) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel-22qDup", "NoASD-ASD", "CNV:ASD")
anova.diagnosis.deldup.asd.noasd$Module <- colnames(ME.genes.DelDup.ASD)
write.xlsx(anova.diagnosis.deldup.asd.noasd, "./celltypesWGCNA.ME.anova.diagnosis.deldup.asd.noasd_012621.xlsx")

#Nominal P vals
anova.diagnosis.deldup.asd.noasd.nogrey.nom.p <- map(ME.genes.DelDup.ASD.nogrey, EigengeneANOVA.deldup.asd.noasd.nogrey, pheno.wME.rmmed.22qDelDup.update.ASD) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.deldup.asd.noasd.nogrey.nom.p) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel-22qDup", "NoASD-ASD", "CNV:ASD")
anova.diagnosis.deldup.asd.noasd.nogrey.nom.p$Module <- colnames(ME.genes.DelDup.ASD.nogrey)
write.xlsx(anova.diagnosis.deldup.asd.noasd.nogrey.nom.p, "./celltypesWGCNA.ME.anova.diagnosis.deldup.asd.noasd.nogrey.nom.p_012221.xlsx")

anova.diagnosis.deldup.asd.noasd.nom.p <- map(ME.genes.DelDup.ASD, EigengeneANOVA.deldup.asd.noasd, pheno.wME.rmmed.22qDelDup.update) %>% 
  reduce(rbind) %>% as_tibble %>%
  mutate_all(signif, digits = 3)
colnames(anova.diagnosis.deldup.asd.noasd.nom.p) <- c("Age.at.Collection.Years", "Sex2", "RIN2", "22qDel-22qDup", "NoASD-ASD", "CNV:ASD")
anova.diagnosis.deldup.asd.noasd.nom.p$Module <- colnames(ME.genes.DelDup.ASD)
write.xlsx(anova.diagnosis.deldup.asd.noasd.nom.p, "./celltypesWGCNA.ME.anova.diagnosis.deldup.asd.noasd.nom.p_012221.xlsx")


### Convert to long form for plotting ###
pheno.wME.rmmed.22qDelDup.update.ASD.long <- gather(pheno.wME.rmmed.22qDelDup.update.ASD, Module, Expression, MEyellow:MEcyan, factor_key = TRUE)
pheno.wME.rmmed.22qDelDup.update.ASD.long$Module <- str_replace_all(pheno.wME.rmmed.22qDelDup.update.ASD.long$Module, "ME", "")
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce <-  pheno.wME.rmmed.22qDelDup.ASD.long[!(pheno.wME.rmmed.22qDelDup.ASD.long$Module == "grey"), ]
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce <- left_join(pheno.wME.rmmed.22qDelDup.ASD.long.reduce, Kangmodule.key)
#pheno.wME.rmmed.22qDelDup.ASD.long.reduce$ModuleNumberName <- factor(pheno.wME.rmmed.22qDelDup.ASD.long.reduce$ModuleNumberName, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18"))

#####  Create the control 0 column and then set values to names in the ASD.CNV.Groups  ####
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0 <- pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0[pheno.wME.rmmed.22qDelDup.update.ASD.long$Diagnosis22q == "control"] <- "0"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups <- "control"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups[pheno.wME.rmmed.22qDelDup.update.ASD.long$Diagnosis22q == "22qDel" & pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0 == "0"] <- "22qDel.NoASD"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups[pheno.wME.rmmed.22qDelDup.update.ASD.long$Diagnosis22q == "22qDel" & pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0 == "1"] <- "22qDel.ASD"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups[pheno.wME.rmmed.22qDelDup.update.ASD.long$Diagnosis22q == "22qDup" & pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0 == "0"] <- "22qDup.NoASD"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups[pheno.wME.rmmed.22qDelDup.update.ASD.long$Diagnosis22q == "22qDup" & pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.DIAGNOS.ANY.TIME.CTL.0 == "1"] <- "22qDup.ASD"
pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups <- factor(pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups, levels = c("control", "22qDel.NoASD", "22qDel.ASD", "22qDup.NoASD", "22qDup.ASD"))

### Keep all except for controls ####
pheno.wME.rmmed.22qDelDup.update.ASD.long.noctl <- pheno.wME.rmmed.22qDelDup.update.ASD.long[pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups == "22qDel.NoASD" | pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups == "22qDel.ASD" | pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups == "22qDup.NoASD" | pheno.wME.rmmed.22qDelDup.update.ASD.long$ASD.CNV.Groups == "22qDup.ASD",]

#### Plot the ME box plots as a facet ###
p <- ggplot(pheno.wME.rmmed.22qDelDup.update.ASD.long.noctl, aes(x = ASD.CNV.Groups, y = Expression, color = ASD.CNV.Groups)) + 
  #intialize ggplot object with data
  #geom_boxplot(width = 0.25, outlier.color=NA, aes(x=ASD.CNV.Groups, y=FNBP1), fill="white") + # y=Num.Genes
  # geom_violin() + #stat = "identity", position="dodge"
  #scale_y_continuous(limits = c(20,135), breaks = c(25, 40, 55, 70, 85, 100, 115, 130)) + #scale_x_continuous(expand = c(0,0))
  facet_wrap(~ Module, ncol = 5) +
  geom_boxplot(width = 0.7) + #Make boxplot
  geom_jitter(width = 0.1) + #Show points
  ylab("Module Eigengene") +
  xlab("Group") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5, color = "black"),
        panel.spacing = unit(2, "lines"),
        #panel.border = element_blank(),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(size = 16, face = "bold", color = "black", angle = 45, hjust=0.6, vjust=0.6), 
        #margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(size=16, color = "black", angle = 0),
        axis.title.x = element_text(size=30, color = "black", angle = 0),
        axis.title.y = element_text(size=30, color = "black", angle = 90),
        legend.position = "none") + #Hide legend
  scale_color_manual(values = c("red", "red4", "green", "green4")) #Override default color scheme
#legend.title=element_text(size=30),
#legend.text=element_text(size=30),
#legend.position = "none")
#legend.margin = margin(30,0,30,30),
#legend.key.width = unit(2.5, "line"),
#legend.key.size = unit(2.8, "lines"))
CairoPNG("ME_plot_CNVASD_facet_012621.png", height = 1400, width = 1400, bg = "transparent")
print(p)
dev.off()
