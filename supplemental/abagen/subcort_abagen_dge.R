# C Schleifer 9/30/22
# this script reads and plots Allen Human Brain Atlas gene expression data previously mapped to CAB-NP atlas surface with abagen
# also reads RSFA values calculated by (which script?) and tests imaging transcriptomic correlations
# cortical parvalbumin and somatostatin expression are visualized as examples
# https://abagen.readthedocs.io/en/stable/user_guide/parcellations.html
# https://github.com/ColeLab/ColeAnticevicNetPartition

####################################################################################################################################################
### set up workspace
####################################################################################################################################################

# clear workspace
rm(list = ls(all.names = TRUE))

# use SSHFS to mount hoffman2 server (download SSHFS for mac: https://osxfuse.github.io/)
# TODO: set hoffman2 username
uname <- "schleife"
# set local path to mount server
hoffman <- "~/Desktop/hoffman_mount"
# create directory if needed 
if(!file.exists(hoffman)){dir.create(hoffman)}
# make string to run as system command
mntcommand <- paste0("umount -f ", hoffman,"; sshfs ",uname,"@hoffman2.idre.ucla.edu:/u/project/cbearden/data ",hoffman)
# if hoffman directory is empty, use system command and sshfs to mount server, if not empty assume already mounted and skip
if(length(list.files(hoffman)) == 0){system(mntcommand)}else{print(paste(hoffman,"is not empty...skipping SSHFS step"))}


# list packages to load
packages <- c("devtools","conflicted","here","magrittr", "dplyr", "tidyr", "ggplot2","ggpubr","RColorBrewer", "ciftiTools","tableone", "data.table", "reshape2","neuroCombat","limma")

# install packages if not yet installed
# note: ciftiTools install fails if R is started without enough memory on cluster (try 16G)
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# use the filter function from dplyr, not stats
conflict_prefer("filter", "dplyr")

# get path to project repo directory
#project <- here()
project <- "/Users/charlie/Dropbox/github/22q_subcort_volumes/"
print(paste("Project directory:", project))



####################################################################################################################################################
### read AHBA data and compute DE between thalamus MDm and LGN
####################################################################################################################################################
# use limma package for differential gene expression: https://bioconductor.org/packages/release/bioc/html/limma.html
# example paper using method: https://www.nature.com/articles/s41380-022-01489-8
# other example paper: https://www.nature.com/articles/s41467-018-03811-x
# https://github.com/HolmesLab/CORTICOSTRIATAL_HOLMES/blob/master/scripts/function_library.R
# https://github.com/HolmesLab/CORTICOSTRIATAL_HOLMES/tree/master


# read individual subject abagen data
abagen_dir <- file.path(project, "abagen")
abagen_names <- list.files(abagen_dir, pattern="thal_native_expression*")
abagen_files <- lapply(file.path(abagen_dir, abagen_names), read.csv, header=T, sep=",", row.names="label")
abagen_ids <- abagen_names %>% gsub("thal_native_expression.","",.) %>% gsub(".csv","",.)
names(abagen_files) <- abagen_ids
# transpose so that rows=genes, columns=brain regions
abagen_t <- lapply(abagen_files, t)

# read abagen info file
thal_info_abagen <- read.csv(file.path(abagen_dir,"thal_info_abagen.csv"), header=T, sep=",")
thal_info_abagen$label_edit <- thal_info_abagen$label %>% gsub("\\(","_",.) %>% gsub("\\)","",.) %>% gsub("-","_",.)

# function to rename columns with region names
col_id_to_name <- function(df, info, donors){
  names <- data.frame(id=colnames(df))
  merged <- merge(names, info[,c("id","label_edit")], by="id")
  colnames(df) <- merged$label_edit
  return(df)
}

abagen_roi_cols <- lapply(abagen_t, function(x) col_id_to_name(df=x, info=thal_info_abagen)) 

# create dataframe with two columns per donor for two regions being compared 
abagen_select_cols <- function(list, cols){
  out <- lapply(list, function(x) x[,cols]) %>% do.call(cbind,.) %>% as.data.frame
  return(out)
}

abagen_mdm_lgn <- abagen_select_cols(list=abagen_roi_cols, cols=c("Left_MDm","Left_LGN"))

# add donor ID to column name
donors <- rep(abagen_ids, times=2)
donors <- donors[order(match(donors, abagen_ids))]
colnames(abagen_mdm_lgn) <- paste0(colnames(abagen_mdm_lgn),"_", donors)

# manually drop columns with NA data 
# TODO: change if using different regions
# rows corresponding to genes and columns to samples
abagen_mdm_lgn_subset <- subset(abagen_mdm_lgn, select=-c(Left_LGN_10021,Left_MDm_15697))

# consider alternative normalization methods
# gil amy gen: "variance stabilized normalization" first normalize variance across genes and then robust line regression
# log transform, first adding 1 to all data to avoid log(0)
exp <- log(abagen_mdm_lgn_subset)
#exp <- abagen_mdm_lgn_subset

# get donor block from column names by removing everything before final underscore
donor_block <- colnames(exp) %>% gsub("^.*\\_","",.) %>% as.numeric

# create design with rows corresponding to samples and columns to coefficients to be estimated.
# first designate MDm as 0 and LGN as 1
sample_grouping <- colnames(exp) %>% gsub("_[^_]+$", "", .) %>% gsub("Left_MDm",0,.) %>% gsub("Left_LGN",1,.) %>% as.numeric
design <- cbind(Grp1=1,Grp2vs1=sample_grouping) 

# DE
dupcor <- duplicateCorrelation(exp,design=design,block=donor_block)
fit <- lmFit(exp, design=design, block=donor_block, correlation=dupcor$consensus)
fit.b <- eBayes(fit)
fit.b$genes <- row.names(exp)

volcanoplot(fit.b,coef=2,highlight=0, names=fit.b$genes)

top_genes <- topTable(fit.b, coef=2, number=10, adjust.method="BH", sort.by = "p")
top_genes

# look into other post-mortem thalamic studies
# check model with statistician
# try different normalization, look into outliers
# look into preproc
# look within donors
# restrict to hypothesis gene set

########################################################################################
# example with sim data from lmFit help
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
#options(digits=3)

# Ordinary fit
fit <- lmFit(y,design)
fit <- eBayes(fit)

# correlated arrays
block <- c(1,1,2,2,3,3)
dupcor <- duplicateCorrelation(y,design,block=block)
dupcor$consensus.correlation
fit3 <- lmFit(y,design,block=block,correlation=dupcor$consensus)
fitE <- eBayes(fit3)
volcanoplot(fitE,coef=2,highlight=2)
volcanoplot(fit,coef=2,highlight=2)
