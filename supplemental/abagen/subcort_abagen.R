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
packages <- c("devtools","conflicted","here","magrittr", "dplyr", "tidyr", "ggplot2","ggpubr","RColorBrewer", "ciftiTools","tableone", "data.table", "reshape2","neuroCombat")

# install packages if not yet installed
# note: ciftiTools install fails if R is started without enough memory on cluster (try 16G)
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# install neuroComBat from github 
# https://github.com/Jfortin1/neuroCombat_Rpackage
#install_github("jfortin1/neuroCombatData")
#install_github("jfortin1/neuroCombat_Rpackage")

# use the filter function from dplyr, not stats
conflict_prefer("filter", "dplyr")

# get path to project repo directory
project <- here()
print(paste("Project directory:", project))

# set up connectome workbench path for ciftiTools
# https://www.humanconnectome.org/software/get-connectome-workbench
# local wbpath (edit this path if workbench is installed in another location, e.g. on hoffman: /u/project/CCN/apps/hcp/current/workbench/bin_rh_linux64/)
# TODO: edit if necessary
wbpath <- "/Applications/workbench/bin_macosx64/"
ciftiTools.setOption("wb_path", wbpath)

# set local path to SSHFS mount point for hoffman:/u/project/cbearden/data
# TODO: edit if necessary
hoffman <- "~/Desktop/hoffman_mount/"

# load rgl for ciftiTools visualization
# may require XQartz v2.8.1 to be installed locally
if(!require('rgl', quietly=TRUE)){install.packages('rgl')}
rgl::setupKnitr()
rgl::rgl.open(); rgl::rgl.close()

####################################################################################################################################################
### read AHBA and CAB-NP data, define plotting functions, and visualize specific genes
####################################################################################################################################################

# load parcellated AHBA data
# this csv is the output of abagen.get_expression_data() python function to extract AHBA expression from CAB-NP surface atlas
# https://abagen.readthedocs.io/en/stable/user_guide/expression.html
ahbaSurfCABNP <- read.csv(file.path(project,"CAB-NP_surface_abagen_expression.csv"), header=T, sep=",")
ahbaSurfCABNP$label2 <- ahbaSurfCABNP$label
# load AHBA extracted from separated CAB-NP volume in MNI space (no cortical ROIs)
ahbaVolCABNP <- read.csv(file.path(project,"CAB-NP_subcort_abagen_expression.csv"), header=T, sep=",")
ahbaVolCABNP$label2 <- ahbaVolCABNP$label

# load CAB-NP network parcellation
# https://github.com/ColeLab/ColeAnticevicNetPartition
ji_key <- read.table(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"),header=T)
ji_net_keys <- ji_key[,c("NETWORKKEY","NETWORK")] %>% distinct %>% arrange(NETWORKKEY)
# read cifti with subcortical structures labeled 
xii_Ji_parcel <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"), brainstructures = "all")
xii_Ji_network <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dscalar.nii"), brainstructures = "all")
# read only surface parcels
xii_Ji_parcel_surf <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"), brainstructures = c("left", "right"))

# function to take xifti atlas (with ROIs denoted by unique values) and return list of xifti matrix indices by brain structure for each ROI
get_roi_atlas_inds <- function(xii){
  # get all unique roi labels
  mxii <- as.matrix(xii)
  vals <- mxii %>% unique %>% sort %>% as.numeric
  # get brain structures to iterate through (cortex_left, cortex_right, subcort)
  xiinames <- names(xii$data)
  # for each ROI, get the indices for each brain structure
  # output is a nested list for each ROI (with "r_" added as prefix) and each brain structure, containing an array of xifti indices corresponding to the ROI in a given brain structure
  out <- lapply(setNames(vals,paste0("r_",vals)), function(v) lapply(setNames(xiinames,xiinames), function(n) which(xii$data[[n]] == v)))
  return(out)
}

# function to create xifti for plotting ROI values on a brain atlas
# input atlas xifti and data frame with at least two cols corresponding to ROI IDs (roi_col) and output values (val_col) e.g. gene expression or functional connectivity
# output modified atlas xifti with ROI IDs replaced with output values (for visualization)
atlas_xifti_new_vals <- function(xii, df, roi_col, val_col){
  # get list of xifti indices for each ROI
  inds <- get_roi_atlas_inds(xii)
  # create blank xii from atlas
  xii_out <- xii
  for (struct in names(xii_out$data)){
    if (!is.null(xii_out$data[[struct]])){
      xii_out$data[[struct]] <- as.matrix(rep(NA,times=nrow(xii_out$data[[struct]])))
    }
  }
  # create new column named roilabel from roi_col
  df$roilabel <- df[,roi_col]
  # for each roi in xii, set all relevant vertices to value from val_col based on roi_col
  for (roi in names(inds)){
    print(roi)
    # get value for roi
    out_val <- as.numeric(filter(df[,c("roilabel",val_col)], roilabel==gsub("r_","",roi))[val_col])
    print(out_val)
    # loop through brain structures, if ROI has any indices in a structure, set those to the output value
    for (struct in names(inds[[roi]])){
      roi_inds <- inds[[roi]][[struct]]
      l <- length(roi_inds)
      if (l > 0){
        xii_out$data[[struct]][roi_inds] <- out_val
      }
    }
  }
  return(xii_out)
}




#########
### get 22q11 gene expression for hip and thal regions
#########

# genes from https://academic.oup.com/cercor/article/31/7/3285/6150031 supplement
forsyth_genes <- c("DGCR8", "AIFM3", "SCARF2", "CLDN5", "DGCR2", "P2RX6", "TANGO2", "RANBP1", "HIRA", "UFD1", "ARVCF", "MED15", "COMT", "PRODH", "SLC25A1", "GNB1L", "GP1BB", "MRPL40", "KLHL22", "RIMBP3", "PI4KA", "RTN4R", "C22orf39", "SEPT5", "DGCR6", "SLC7A4", "DGCR6L", "SNAP29")

# n=3 not in AHBA volume CABNP
forsyth_genes_missing <- forsyth_genes[which(!forsyth_genes %in% names(ahbaVolCABNP))]
forsyth_genes_use <- forsyth_genes[which(forsyth_genes %in% names(ahbaVolCABNP))]

# get main structure label from xifti metadata for each subcortical roi (361:718)
#sc_parc_labs <- lapply(361:718, function(r) xii_Ji_parcel$meta$subcort$labels[which(xii_Ji_parcel$data$subcort==r)])
#sc_parc_labs <- lapply(sc_parc_labs, droplevels)
#sc_parc_labs_char <- data.frame(INDEX=361:718,structure=do.call(rbind,lapply(sc_parc_labs, levels)))
#ji_key_sc <- merge(x=filter(ji_key, INDEX >= 361 & INDEX <= 718), y=sc_parc_labs_char, by="INDEX")

# get key for subcortex only
ji_key_sc <- filter(ji_key, INDEX >= 361 & INDEX <= 718)

# get structure name from LABEL column
ji_key_sc$structure <- do.call(rbind,strsplit(ji_key_sc$LABEL, split="_"))[,2]

# get all right hippocampal parcels
ji_key_sc_rhipp <- filter(ji_key_sc, structure=="R-Hippocampus")

# get all right thal frontoparietal parcels
ji_key_sc_rthal_fpn <- filter(ji_key_sc, structure=="R-Thalamus" & NETWORK =="Frontoparietal")
# single parcel for right thal fpn
ind_rthal_fpn <- ji_key_sc_rthal_fpn$INDEX
# 22q gene expression for right thal fpn
exp22q_rthal_fpn <- filter(ahbaVolCABNP, label==ind_rthal_fpn)[,forsyth_genes_use]
exp22q_rthal_fpn <- apply(exp22q_rthal_fpn,2,mean)
# 22q gene expression for right hipp
exp22q_rhipp <- filter(ahbaVolCABNP, label %in% ji_key_sc_rhipp$INDEX, !is.na(DGCR8))[,forsyth_genes_use]
exp22q_rhipp_mean <- apply(exp22q_rhipp,2,mean)

# make data frame for plotting
rhipp <- data.frame(exp=exp22q_rhipp_mean, gene=names(exp22q_rhipp_mean), structure="hipp_mean")
rthal <- data.frame(exp=exp22q_rthal_fpn, gene=names(exp22q_rthal_fpn), structure="thal_fpn")
thal_hip <- rbind(rthal, rhipp)

ggplot(thal_hip, aes(exp, fill=structure))+
  geom_density(kernel="gaussian", alpha=0.3)+
  scale_fill_manual(values=c("lightblue","red"))+
  geom_vline(xintercept = mean(exp22q_rhipp_mean), lty="dashed")+
  geom_vline(xintercept = mean(exp22q_rthal_fpn), lty="dashed")+
  theme_classic()

# use limma package for differential gene expression: https://bioconductor.org/packages/release/bioc/html/limma.html
# example paper using method: https://www.nature.com/articles/s41380-022-01489-8
