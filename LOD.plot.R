library(data.table)
library(tidyverse)
library(cowplot)

'%ni%' <- function(x,y)!('%in%'(x,y))

##set plotting variables
pd <- position_dodge(0.5)
fontSize = 18
pointSize = 6

##take in files
#setwd('/Volumes/')
dilution_table <- fread('cgi/scratch/fbeaudry/msi_test/LOD/LOD.msi')
names(dilution_table) <- c("sample","dilutionFactor","q0","q1","median","q3","q4")

##edit sample dilution factor
#multiply dilution by coverage and then by purity, for dilution factors with one sig.fig. divide by 10, for D.F. of 2 S.F. divide by 100
dilution_table$dilutionFactor_corr[dilution_table$sample == "OCT_010472" ] <- (dilution_table$dilutionFactor[dilution_table$sample == "OCT_010472" ]/10) * 103 * 0.71
dilution_table$dilutionFactor_corr[dilution_table$sample == "OCT_010472" & dilution_table$dilutionFactor %in% c(25)] <- (dilution_table$dilutionFactor[msi$sample == "OCT_010472"  & dilution_table$dilutionFactor %in% c(25)]/100) * 103 * 0.71

dilution_table$dilutionFactor_corr[dilution_table$sample == "OCT_010676" ] <- 
  (dilution_table$dilutionFactor[dilution_table$sample == "OCT_010676" ]/10) * 88 * 0.62
dilution_table$dilutionFactor_corr[dilution_table$sample == "OCT_010676" & dilution_table$dilutionFactor %in% c(35,75)] <- 
  (dilution_table$dilutionFactor[dilution_table$sample == "OCT_010676" & dilution_table$dilutionFactor  %in% c(35,75)]/100 ) * 88 * 0.62

##plot
ggplot(dilution_table,
       aes(x=dilutionFactor_corr,y=median,color=sample )) + 
  geom_hline(yintercept = 3.5,color="grey",linetype="dotted")+
  geom_hline(yintercept = 5,color="grey",linetype="dashed")+
  
  geom_vline(xintercept = 24,color="grey")+ #24X coverage is our minimum _tumor_ coverage for either assay
  
  geom_point( shape=1, size=pointSize) + 
  geom_line( ) + 
  
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.1) +
  theme_bw(base_size = fontSize) + 
  labs(x="Post-dilution Tumor Coverage (Purity * Coverage)",y="Score",color="Sample") + 

  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
