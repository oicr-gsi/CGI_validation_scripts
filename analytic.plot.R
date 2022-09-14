#library(data.table)
library(tidyverse)
library(cowplot)

'%ni%' <- function(x,y)!('%in%'(x,y))

make_cancerTypeMatrix <- function(sample_Ctype){
  #downcast matrix into cancer tissue columns
  sample_Ctype_mat <- reshape2::dcast(sample_Ctype,formula=Sample+Project+Study~Cancertissue,value.var ="Cancertissue" ) 
  
  #remove the sample info columns for binarization
  sample_Ctype_mat_01 <- sample_Ctype_mat[ , !(names(sample_Ctype_mat) %in% c("Sample","Project","Study"))]
  
  
  sample_Ctype_mat_01[!is.na(sample_Ctype_mat_01)] <- 1
  sample_Ctype_mat_01[is.na(sample_Ctype_mat_01)] <- 0
  
  #add back sample info
  sample_Ctype_mat_01_lab <- cbind.data.frame(sample_Ctype_mat[,c(1:3)],sample_Ctype_mat_01)
  sample_Ctype_mat_01_melt <- reshape2::melt(sample_Ctype_mat_01_lab,id.vars=c("Sample","Project","Study"))
  
  #only keep positive rows
  sample_Ctype_mat_01_melt <- sample_Ctype_mat_01_melt[sample_Ctype_mat_01_melt$value == 1,]
  return(sample_Ctype_mat_01_melt)
}

#set file names
cohort_file <- '~/Documents/data/Validation_cohort_WGS.txt'
mutations_file <-   '~/Documents/data/cohort_zero.mutations.txt'
msi_file <- '/Volumes/cgi/scratch/fbeaudry/msi_test/all.msi'

inconclusive_cutoff = 3.5
MSI_cutoff = 5

pd <- position_dodge(0.5)
fontSize = 18
pointSize = 6

#import files to dataframes#
cohort <- read.table(cohort_file,header = T)
mutations <- read.table(mutations_file,header = T)
msi <- read.table(msi_file,header = T)

#cancer type matrix
sample_Ctype <- cohort %>% select(Sample,Cancertissue,Project,Study)
sample_Ctype_mat_01_melt <- make_cancerTypeMatrix(sample_Ctype)

##process MSI dataframe for consistancy
msi <- separate(data=msi,col=sample,into=c("Study","Sample"),sep ="_",remove=FALSE,extra="drop")

#add MSI call from score (total VS median??)
msi$msiCall <- "MSS"
msi$msiCall[msi$total >= inconclusive_cutoff] <- "Inconclusive"
msi$msiCall[msi$total >= MSI_cutoff] <- "MSI"

####MMR gene status####
MMRD_vars <- mutations %>% 
  filter(gene %in% c("MLH1", "MSH2", "MSH6")) %>% 
  merge(sample_Ctype_mat_01_melt,by.y="Sample",by.x="donor",all=TRUE)

#make dummy variable for mutation panel in plot 
MMRD_vars$gene[is.na(MMRD_vars$gene)] <- "MSH6"

MMRD_vars$mutationCategory[MMRD_vars$mutationCategory %in% c("DEL","INS","SNP")] <- "small"
MMRD_vars$MutationCategory <-  factor(MMRD_vars$mutationCategory, levels = c("small","LOH","CNV"))
MMRD_vars_smallIDs <- MMRD_vars$Sample[MMRD_vars$MutationCategory == "small"]

MMRD_vars$mutationType[MMRD_vars$mutationType %in% c("LOH","Deletion")] <- NA

MMRD_vars$MutationSymbol[MMRD_vars$MutationCategory == "small"] <- NA
MMRD_vars$MutationSymbol[MMRD_vars$MutationCategory == "CNV"] <- "DEL"
MMRD_vars$MutationSymbol[MMRD_vars$MutationCategory == "LOH"] <- "LOH"
MMRD_vars$MutationSymbol[MMRD_vars$Sample %ni% MMRD_vars_smallIDs & MMRD_vars$MutationSymbol == "LOH"] <- NA

####plot####

##sample score as factor order/level
sample_order <- msi$sample[  order(  msi$total)]
msi$Sample                      <-  factor(msi$sample,  levels = sample_order)
sample_Ctype_mat_01_melt$sample <- factor(sample_Ctype_mat_01_melt$Sample, levels = sample_order)
MMRD_vars$sample                <-  factor(MMRD_vars$donor, levels = sample_order)
cohort$sample                   <-  factor(cohort$Sample, levels = sample_order)

MSI_plot <- 
  ggplot(msi,aes(x=Sample,color=msiCall)) + 
    geom_hline(yintercept = 3.5,color="grey",linetype="dotted")+
    geom_hline(yintercept = 5,color="grey",linetype="dashed")+
    
    geom_point(aes(y=total), shape=1, size=pointSize) + 
    geom_point(aes(y=q50), shape=4, size=pointSize) + 
    
    geom_errorbar(aes(ymin=q25, ymax=q75), width=.1) +
    
    theme_bw(base_size = fontSize) + labs(x="Sample",y="MSI Score (% of unstable MS)",color="MSI call") + 
    theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))  +
  
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank()
          )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

MMR_info_plot <- 
  ggplot(MMRD_vars %>% filter(sample %in% msi$sample) %>% arrange(MutationCategory), aes(x=sample,y=gene,fill=mutationType)) + 
      
    geom_tile(aes(fill=mutationType)) +
    geom_point(aes(size=MutationSymbol),size=pointSize,color="white")+

    labs(x="Sample",y="Gene",fill="Mutation") + 
  
    scale_shape_manual(values=c(13,1))   +
    scale_fill_discrete(na.value=NA,breaks=c("Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Splice_Site")) +
  
    theme_bw(base_size = fontSize) + 
    theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) 
  
cancer_type_plot <- 
  ggplot(sample_Ctype_mat_01_melt  %>% filter(sample %in% msi$sample), aes(x=sample, y=variable, fill= Study)) + 
    
    geom_tile() + 
    
    labs(fill="Study",y="Tissue",x="Sample") + 
    guides(alpha="none") +
  
    theme_bw(base_size = fontSize) + 
    theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
    theme(legend.background = element_rect(colour = NA,fill=NA))

TMB_plot <- 
  ggplot(cohort, aes(x=sample, y="TMB", fill= TMB)) + 
    geom_tile() + 
  
    labs(fill="",y=" ") +
  
    scale_fill_gradient2(
      high = "red",
      mid = "white",
      low = "blue",
      midpoint = 3
    ) +
  
    theme_bw(base_size = fontSize) + 
    theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks=element_blank()) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) 

now <- format(Sys.time(), "%d%b%Y")

png(paste0("~/Documents/data/msi/msi.",now,".png"), width = 1200, height = 1000)

plot_grid(MSI_plot, MMR_info_plot, TMB_plot,cancer_type_plot, ncol = 1, rel_heights = c(2,  1,0.5,1.5),align = 'v')

dev.off()
 
####clustering analysis####

library(mclust)

density_model <- densityMclust(msi$q50)
summary(density_model)
m_model <- Mclust(msi$q50,G=2)
summary(m_model)

msi$cluster_class <- m_model$classification


