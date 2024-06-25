
#FUNCTION WITH ALL THE DIFFERENT CONDITIONS/TREATMENT/SEX/GENOTYPE AND SAMPLES NAMES

Sample_Conditions <- function(RawProteinDataWithSampleNames){
  
  targets <- RawProteinDataWithSampleNames %>% 
    dplyr::select(starts_with('F'),starts_with('M'))
  
  targets <- data.frame(sample = colnames(targets)) %>% 
    mutate(Sex = case_when(startsWith(sample,"F") ~ "Female",
                           startsWith(sample,"M") ~ "Male"
    ))
  targets <- targets %>%
    mutate(Condition = case_when(
      startsWith(sample, "F WT KCl") ~ "KCl",
      startsWith(sample, "M WT KCl") ~ "KCl",
      startsWith(sample, "F KO KCl") ~ "KCl",
      startsWith(sample, "M KO KCl") ~ "KCl",
      startsWith(sample, "F WT IS") ~ "IS",
      startsWith(sample, "M WT IS") ~ "IS",
      startsWith(sample, "F KO IS") ~ "IS",
      startsWith(sample, "M KO IS") ~ "IS"
    ))
  targets <- targets %>%
    mutate(Pool = case_when(
      startsWith(sample, "F WT KCl") ~ "F WT",
      startsWith(sample, "M WT KCl") ~ "M WT",
      startsWith(sample, "F KO KCl") ~ "F KO",
      startsWith(sample, "M KO KCl") ~ "M KO",
      startsWith(sample, "F WT IS") ~ "F WT",
      startsWith(sample, "M WT IS") ~ "M WT",
      startsWith(sample, "F KO IS") ~ "F KO",
      startsWith(sample, "M KO IS") ~ "M KO"
    ))
  
  targets <- targets %>%
    mutate(Genotype = case_when(
      startsWith(sample,"F WT") ~ "WT",
      startsWith(sample,"M WT") ~ "WT",
      startsWith(sample,"F KO") ~ "KO",
      startsWith(sample,"M KO") ~ "KO",
    ))
  targets <- targets %>%
    mutate(Treatment_groups = case_when(
      startsWith(sample, "F WT KCl") ~ "F WT KCl",
      startsWith(sample, "M WT KCl") ~ "M WT KCl",
      startsWith(sample, "F KO KCl") ~ "F KO KCl",
      startsWith(sample, "M KO KCl") ~ "M KO KCl",
      startsWith(sample, "F WT IS") ~  "F WT IS",
      startsWith(sample, "M WT IS") ~ "M WT IS",
      startsWith(sample, "F KO IS") ~ "F KO IS",
      startsWith(sample, "M KO IS") ~ "M KO IS"
    ))
  return(as.data.frame(targets))
}

#DATA quality check. We perform basic QC plots in the data to see if the data is OK to proceed to the next steps
QualityCheck <- function(DataToCheck){
  
  colsums <- colSums(log2(DataToCheck), na.rm = TRUE)
par(mfrow = c(2,2))
barplot(colsums)
#sd variations across columns
sd_by_col<- log2(DataToCheck)%>% 
  summarise_if(is.numeric, sd, na.rm = TRUE)
hist(as.matrix(sd_by_col))
#sd_variations across rows
sd_by_row <- as.data.frame(rowSds(as.matrix(log2(DataToCheck))))
names(sd_by_row) <- c("sd_by_prot")
hist(as.matrix(sd_by_row))

#PLOT DENSITIES 
plotDensities(DataToCheck, legend=FALSE) # When we look at density plots of the intensities, we see that most intensities are low and that there is a long tail to the right. It is more useful to explore the data upon log transformation.
#plotDensities(log2(DataToCheck), legend=FALSE) # When we look at the densities we see small shifts in location and shape of the density distributions. So normalisation will be required
  
}

#DATA SUBSETTING

SubsetLogCountsNoNAsDF <- function (Df){
    Females <- Df %>% 
    dplyr::select(starts_with('F'))%>% 
    log2()%>%
    drop_na()
  
    Males <- Df %>% 
      dplyr::select(starts_with('M'))%>% 
      log2()%>%
      drop_na()
    
    BothSexes <- Df %>% 
      dplyr::select(starts_with('F'), starts_with('M'))%>% 
      log2()%>%
      drop_na()
    
   # ListOfSexSubsets <- tibble(groups = c("Females", "Males", "BothSexes"),
   #                             subsets = list(
    #                           tibble(Females),
     #  #                        tibble(Males), 
         #                      tibble(BothSexes)))

   # ListOfSexSubsets <- list(Females, Males , BothSexes)
    
   # names(ListOfSexSubsets) <- c("Females", "Males" , "BothSexes")
    
    return(Females, Males, BothSexes)
}

#VIOLIN PLOTS

PlotViolin <- function(Df){
Df <- Df %>% as.data.frame() %>% rownames_to_column("id")
Df_violin <- melt.data.frame(Df,"id")

Df_violin2 <- as.data.frame(as.character(Df_violin$variable))
names(Df_violin2) <- "variable"
Df_violin2 <- Df_violin2 %>% 
  mutate(Condition = case_when(
    startsWith(variable, "F WT KCl") ~ "KCl",
    startsWith(variable, "F KO KCl") ~ "KCl",
    startsWith(variable, "F WT IS") ~ "IS",
    startsWith(variable, "F KO IS") ~ "IS",
    startsWith(variable, "M WT KCl") ~ "KCl",
    startsWith(variable, "M KO KCl") ~ "KCl",
    startsWith(variable, "M WT IS") ~ "IS",
    startsWith(variable, "M KO IS") ~ "IS" 
  ),
  Genotype = case_when(
    startsWith(variable,"F WT") ~ "WT",
    startsWith(variable,"F KO") ~ "KO",
    startsWith(variable,"M WT") ~ "WT",
    startsWith(variable,"M KO") ~ "KO"
  ),
  Treatment_groups = case_when(
    startsWith(variable, "F WT KCl") ~ "F WT KCl",
    startsWith(variable, "F KO KCl") ~ "F KO KCl",
    startsWith(variable, "F WT IS") ~  "F WT IS",
    startsWith(variable, "F KO IS") ~ "F KO IS",
    startsWith(variable, "M WT KCl") ~ "M WT KCl",
    startsWith(variable, "M KO KCl") ~ "M KO KCl",
    startsWith(variable, "M WT IS") ~  "M WT IS",
    startsWith(variable, "M KO IS") ~ "M KO IS"
  )) 

Df_violin <- cbind(Df_violin  , Df_violin2[!names(Df_violin2) %in% "variable"]) 

VIOLIN <- ggplot(Df_violin , aes(x=variable , y=value, fill=Treatment_groups)) + 
  geom_violin(linewidth = 0.25) + 
  scale_fill_manual(values = c( "khaki1", "lightskyblue","gold","dodgerblue4")) + 
  theme(axis.text.y = element_text(size = 25), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_text(size  = 30),legend.text = element_text(size = 15), axis.title.x = element_blank()) 

 return(VIOLIN)   
}

#PCA FUNCTION TO PLOT ALL 3 PCAs FOR THE 3 FIRST PCs + THEIR DENSITY PLOT
PlotPCA_batch <- function(PCA_DF, Results_PCA){
  PC1_PC2 <- ggplot(PCA_DF , aes(PC1,PC2, color=Batch)) + labs(x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC2  (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
    geom_point()+ 
    theme_minimal() + 
    geom_text_repel(aes(label = sample)) +
    #scale_color_manual(values = c("plum2","blue","orange", "green")) +
    theme(legend.position="none")
  PC1_PC3 <-ggplot(PCA_DF , aes(PC1,PC3, color=Batch)) + labs( x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC3 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
    geom_point() + 
    theme_minimal() + 
    geom_text_repel(aes(label = sample)) +
    #scale_color_manual(values = c("plum2","blue","orange","green")) +
    theme(legend.position="none")
  PC3_PC2 <- ggplot(PCA_DF , aes(PC2,PC3, color=Batch)) +
    geom_point()  +
    theme_minimal() +
    geom_text_repel(aes(label = sample)) +
    labs( x = paste("PC2 (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC3 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)")) 
  #scale_color_manual(values = c("plum2","blue","orange", "green"))
  # + stat_ellipse(geom = "polygon", aes(fill = Treatment_groups), alpha = 0.25) (if I want to add ellipses)
  
  #PC COMPONENTS % OF EXPLAINED VARIANCE
  explained <- (Results_PCA$sdev)^2 / sum(Results_PCA$sdev^2)
  explained <- as.data.frame(explained)
  explained$Components <- as.numeric(row.names(explained))
  names(explained)[1] <- "% of explained variance"
  
  boulder<- ggplot(explained, aes(x = Components, y = `% of explained variance` , fill = "orange")) +
    geom_col() +
    theme_minimal() +
    scale_y_continuous(position = "left") +
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_x_continuous(position = "top", breaks = c(1:length(explained[,1])))
  
  #PCA DENSITIES BY COMPONENT 
  
  DENSITY_PC1<- ggplot(PCA_DF, aes(x=PC1)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC1",x="Intensities", y = "Density") 
  
  DENSITY_PC2<- ggplot(PCA_DF, aes(x=PC2)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC2",x="Intensities", y = "Density") 
  
  DENSITY_PC3 <- ggplot(PCA_DF, aes(x=PC3)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC3",x="Intensities", y = "Density") 
  
  
  # Draw Only legend 
  legend_extracted <- cowplot::get_legend(PC3_PC2 )              
  
  # Create new plot window
  grid.newpage()                              
  
  legend <- grid.draw(legend_extracted) 
  
  #plotting everything together
  
  return(ggarrange(DENSITY_PC1, legend, boulder,PC1_PC2, DENSITY_PC2,PC1_PC3, PC3_PC2, DENSITY_PC3, nrow= 3, ncol = 3))
  
}
#PCA FUNCTION TO PLOT ALL 3 PCAs FOR THE 3 FIRST PCs + THEIR DENSITY PLOT
PlotPCA <- function(PCA_DF, Results_PCA){
PC1_PC2 <- ggplot(PCA_DF , aes(PC1,PC2, color=Genotype)) + labs(x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC2  (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
  geom_point()+ 
  theme_minimal() + 
  geom_text_repel(aes(label = sample)) +
  scale_color_manual(values = c("blue","orange")) +
  theme(legend.position="none")
PC1_PC3 <-ggplot(PCA_DF , aes(PC1,PC3, color=Genotype)) + labs( x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC3 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
  geom_point() + 
  theme_minimal() + 
  geom_text_repel(aes(label = sample)) +
  scale_color_manual(values = c("blue","orange")) +
  theme(legend.position="none")
PC3_PC2 <- ggplot(PCA_DF , aes(PC2,PC3, color=Genotype)) +
  geom_point()  +
  theme_minimal() +
  geom_text_repel(aes(label = sample)) +
  labs( x = paste("PC2 (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC3 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
  scale_color_manual(values = c("blue","orange")
  + stat_ellipse(geom = "polygon", aes(fill = Treatment_groups), alpha = 0.25))
# (if I want to add ellipses)

#PC COMPONENTS % OF EXPLAINED VARIANCE
explained <- (Results_PCA$sdev)^2 / sum(Results_PCA$sdev^2)
explained <- as.data.frame(explained)
explained$Components <- as.numeric(row.names(explained))
names(explained)[1] <- "% of explained variance"

boulder<- ggplot(explained, aes(x = Components, y = `% of explained variance` , fill = "orange")) +
  geom_col() +
  theme_minimal() +
  scale_y_continuous(position = "left") +
  theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_continuous(position = "top", breaks = c(1:length(explained[,1])))

#PCA DENSITIES BY COMPONENT 

DENSITY_PC1<- ggplot(PCA_DF, aes(x=PC1)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
  labs(title="PC1",x="Intensities", y = "Density") 

DENSITY_PC2<- ggplot(PCA_DF, aes(x=PC2)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
  labs(title="PC2",x="Intensities", y = "Density") 

DENSITY_PC3 <- ggplot(PCA_DF, aes(x=PC3)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
  labs(title="PC3",x="Intensities", y = "Density") 


# Draw Only legend 
legend_extracted <- cowplot::get_legend(PC3_PC2 )              

# Create new plot window
grid.newpage()                              

legend <- grid.draw(legend_extracted) 

#plotting everything together

return(ggarrange(DENSITY_PC1, legend, boulder,PC1_PC2, DENSITY_PC2,PC1_PC3, PC3_PC2, DENSITY_PC3, nrow= 3, ncol = 3))
 
}

#PLOT BOTH SEXES
PlotPCA_allsexes <- function(PCA_DF, Results_PCA){
  PC1_PC2 <- ggplot(PCA_DF , aes(PC1,PC2, color=Sex)) + labs(x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC2  (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
    geom_point()+ 
    theme_minimal() + 
    scale_color_manual(values = c("blue","orange")) + 
    geom_text_repel(aes(label = sample)) +
    theme(legend.position="none")
  PC1_PC3 <-ggplot(PCA_DF , aes(PC1,PC2, color=Sex)) + 
    labs( title="PCA", 
          x = paste("PC1 (",as.character(round(Results_PCA$sdev[1]^2/sum(Results_PCA$sdev^2)*100)),"%)"), 
          y = paste("PC2 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)")) +
    geom_point() + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12)) + 
    scale_color_manual(values = c("hotpink","navy")) + 
    geom_text_repel(aes(label = sample)) + 
    stat_ellipse(geom = "polygon", aes(fill = Sex), alpha = 0.20)
  PC3_PC2 <- ggplot(PCA_DF , aes(PC2,PC3, color=Sex)) +
    geom_point()  +
    theme_minimal() +
    geom_text_repel(aes(label = sample)) +
    labs( x = paste("PC2 (",as.character(round(Results_PCA$sdev[2]^2/sum(Results_PCA$sdev^2)*100)),"%)"), y = paste("PC3 (",as.character(round(Results_PCA$sdev[3]^2/sum(Results_PCA$sdev^2)*100)),"%)") + 
    stat_ellipse(geom = "polygon", alpha = 0.2))
  # + stat_ellipse(geom = "polygon", aes(fill = Sex), alpha = 0.2))
  
  #PC COMPONENTS % OF EXPLAINED VARIANCE
  explained <- (Results_PCA$sdev)^2 / sum(Results_PCA$sdev^2)
  explained <- as.data.frame(explained)
  explained$Components <- as.numeric(row.names(explained))
  names(explained)[1] <- "% of explained variance"
  
  boulder<- ggplot(explained, aes(x = Components, y = `% of explained variance` , fill = "orange")) +
    geom_col() +
    theme_minimal() +
    scale_y_continuous(position = "left") +
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_x_continuous(position = "top", breaks = c(1:length(explained[,1])))
  
  #PCA DENSITIES BY COMPONENT 
  
  DENSITY_PC1<- ggplot(PCA_DF, aes(x=PC1)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC1",x="Intensities", y = "Density") 
  
  DENSITY_PC2<- ggplot(PCA_DF, aes(x=PC2)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC2",x="Intensities", y = "Density") 
  
  DENSITY_PC3 <- ggplot(PCA_DF, aes(x=PC3)) + geom_density(size=1)+theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5,size = 20),axis.text.y = element_text(size = 12), axis.text.x  = element_text(size = 12), axis.title = element_text(size  = 15)) +
    labs(title="PC3",x="Intensities", y = "Density") 
  
  
  # Draw Only legend 
  legend_extracted <- cowplot::get_legend(PC1_PC3)              
  
  # Create new plot window
  grid.newpage()                              
  
  legend <- grid.draw(legend_extracted) 
  #plotting everything together
  
  # return(ggarrange(DENSITY_PC1, legend, boulder,PC1_PC2, DENSITY_PC2,PC1_PC3, PC3_PC2, DENSITY_PC3, nrow= 3, ncol = 3))
  
  return(PC1_PC3)
}


