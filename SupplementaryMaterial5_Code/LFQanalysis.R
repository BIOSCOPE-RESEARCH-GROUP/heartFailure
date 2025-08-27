#......................          README       ..............................####
# Title: LFQ Proteomic Data Analysis with Feature Selection and Visualization
# Author: Julia García Currás
# Description: This script loads an LFQ proteomic data matrix, computes the mean 
#             between pairs of replicates for each protein, and adjusts the data  
#             structure. Then, five different feature selection (FS) algorithms  
#             are applied to identify the top 20 most relevant proteins. The 
#             script also includes code to visualize the combined results of  
#             these algorithms using a Venn diagram and a clustered heatmap. 
#             Individual predictors and classfication metrics for the most shared 
#             proteins can be found at the end of the script. 
# Input: "in/Data All Proteins.txt" (MaxQuant output)
# Output: "out/featureSelectionResults_LFQ.xlsx", "out/vennDiagram_LFQ.png", "out/heatmap_LFQ.png"
# Dependencies: FSelector, mRMRE, igraph, survival, caret, lattice, ggplot2, 
#               Boruta, dplyr. R version: 4.4.1. 




#...........................................................................####
# Setup ####
rm(list=ls())
graphics.off()
# IMPORTANT! SET YOUR WORKING DIRECTORY TO THIS FOLDER
setwd("your/current/directory")
library(dplyr)
library(Boruta)
library(caret)


#...........................................................................####
# 1. Load data ####
## LFQ data ####
  # Open data from Supplementary Excel
data0 <- read.table(file = "in/Data All Proteins.txt", header = T, sep = "\t")
data0 <- as.data.frame(data0)
  # Select interesting columns
data <- data0 %>% 
  dplyr::select(LFQ_Ctrl_U57_GC1_01_2070:LFQ_HF.rEF_U53_GB5_01_2045, Gene.names) %>%
  dplyr::rename(geneNames = Gene.names)
  # Check duplicated genes names
table(duplicated(data$geneNames))
  # Remove duplicates
data <- data[!duplicated(data$geneNames),]
  # Set rownames and remove the corresponding column
rownames(data) <- data$geneNames
data$geneNames <- NULL


## Design matrix ####
grupos <- unname(sapply(sapply(colnames(data), strsplit, split = "_"), "[[", 2))
individuos <- unname(sapply(sapply(colnames(data), strsplit, split = "_"), "[[", 3))
  # Create data frame with group, sample name and replicate correspondece
designMatrix <- data.frame(
  ID = colnames(data),
  Grupos = factor(grupos, levels = c("Ctrl", "HF.pEF", "HF.rEF")),
  Samples = individuos,
  Replicates = rep(1:2, length(grupos)/2))
  # Design matrix for unique samples
dfGrupos <- designMatrix %>% filter(Replicates == 1) %>% dplyr::select(Samples, Grupos)
colnames(dfGrupos) <- c("Samples", "Groups")
dfGrupos <- dfGrupos[!duplicated(dfGrupos$Samples),]


## Data ####
# Function to average replicates by sample and protein
removeDiscordancesReplicates <- function(dfQ, nReps = 2){
  posFinal <- nReps
  dfResult <- NULL
  for (i in seq(1, ncol(dfQ)-1, nReps)){
    if (i != 1){
      posFinal <- posFinal + nReps
    }  
    dfAux <- dfQ[, i:posFinal] # Selecting replicate rows 
    
    resultProt <- apply(dfAux, 1, function(fila){ # iteration by proteins
      if (any(is.na(fila))){ # NA values not allowed
        return(NA)
      } else {
        return(mean(fila, na.rm = T)) # Average of replicates quantity
      }
    })
    
    if (is.null(dfResult)){ # binding columns
      dfResult <- data.frame(resultProt)
    } else {
      dfResult <- cbind(dfResult, resultProt)
    }
  }
  colnames(dfResult) <- colnames(dfQ)[(seq(nReps, ncol(dfQ), nReps))]
  return(dfResult)
}
  # Apply average to replicates by protein and sample
data <- removeDiscordancesReplicates(dfQ = data)
colnames(data) <- designMatrix[which(designMatrix$ID %in% colnames(data)), "Samples"]




#...........................................................................####
# 2. Adjust data structure ####
df <- data %>% na.omit # Just a checkpoint, Boruta algorithm does not deal with missing data
df <- as.data.frame(t(df))
# Adding group information
df$Grupos <- as.character(dfGrupos[which(dfGrupos$Samples %in% rownames(df)), "Groups"])
table(df$Grupos)
# Removing control group
df <- df %>% 
  filter(Grupos != "Ctrl") %>% 
  dplyr::select(Grupos, everything())
df$Grupos <- factor(df$Grupos, levels = c("HF.pEF", "HF.rEF"))




#...........................................................................####
# 3. Boruta algorithm ####

## Boruta execution ####
set.seed(1234)
bor_en <- Boruta(Grupos ~., data=df, doTrace=2, maxRuns= 500)
bor_en_signif <- names(bor_en$finalDecision[bor_en$finalDecision %in% c("Confirmed")])
# graphical representation of the output
plot(bor_en, xlab="", main="Variable Importance",  las = 2, cex.axis = 0.8, 
     colCode =  c("darkolivegreen4", "orange", "tomato3", "cadetblue4"))
# final test for fixed tentative attributes
Tentative.boruta_en <- TentativeRoughFix(bor_en)
boruta_en.df <- attStats(Tentative.boruta_en)
# Data structuration
boruta_en.df <- boruta_en.df[, c(1,6)]
boruta_en_sort <- boruta_en.df[order(boruta_en.df$meanImp, decreasing = TRUE), ]
boruta_en_sort$variables <- rownames(boruta_en_sort) 
rownames(boruta_en_sort) <- 1:nrow(boruta_en_sort)
# Selecting relevant information to combine results
borutaResults <- boruta_en_sort %>% 
  dplyr::select(variables, meanImp, decision)
colnames(borutaResults) <- c("Proteins", "Importance", "Decision")
# Saving
xlsx::write.xlsx(x = borutaResults, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Boruta", 
                 col.names = T, row.names = F, showNA = F, 
                 append = F)
borutaResults <- borutaResults[, -3]





#...........................................................................####
# 4. Other feature selection algorithms ####

## Data: same as in Boruta algorithm ####
dfFS <- df
g1 = "HF.pEF"
g2 = "HF.rEF"
dfFS <- dfFS %>% dplyr::rename(Groups = Grupos)
dfFS$Groups <- ifelse(dfFS$Groups == g2, 1, 0)
dfFS$Groups <- factor(dfFS$Groups, ordered = T)


## 4.1 Filter methods ####
### 4.1.1 Minimum redundancy maximum relevance ####
library(mRMRe)
set.seed(1234)
dd <- mRMR.data(data = data.frame(dfFS))
seleccion <- mRMR.classic(data = dd, target_indices = c(1), 
                          feature_count = ncol(dfFS)-1)
# Selecting relevant information to combine results
protsFilter <- colnames(dfFS)[as.vector(solutions(seleccion)[[1]])]
scoresProtsFilter <- as.vector(scores(seleccion)[[1]])
names(scoresProtsFilter) <- protsFilter
mrmrResults <- data.frame(Proteins = protsFilter, 
                            Score = scoresProtsFilter)
mrmrResults <- na.omit(mrmrResults)
colnames(mrmrResults) <- c("Proteins", "Importance")
# Saving
xlsx::write.xlsx(x = mrmrResults, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "mRMR", 
                 col.names = T, row.names = F, showNA = F, 
                 append = T)


### 4.1.2 information gain ####
library(FSelector)
set.seed(1234)
info <- FSelector::information.gain(Groups ~ ., data = dfFS)
info$variable <- rownames(info)
# Selecting relevant information to combine results
gainResuls <- info %>% 
  dplyr::select(variable, attr_importance) %>%
  arrange(desc(attr_importance)) %>%
  rename(importance = attr_importance)
colnames(gainResuls) <- c("Proteins", "Importance")
# Saving
xlsx::write.xlsx(x = gainResuls, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Information Gain", 
                 col.names = T, row.names = F, showNA = F, 
                 append = T)


### 4.1.3 Importance measure by ROC (for classification) ####
set.seed(1234)
roc_imp <- caret::filterVarImp(x = dfFS[,2:(ncol(dfFS))], y = dfFS[,1])
# sort the score in decreasing order
roc_imp <- data.frame(cbind(variable = rownames(roc_imp), score = roc_imp[,1]))
roc_imp$score <- as.double(roc_imp$score)
rocResults <- roc_imp[order(roc_imp$score,decreasing = TRUE),]
colnames(rocResults) <- c("Proteins", "Importance")
# Saving 
xlsx::write.xlsx(x = rocResults, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "ROC", 
                 col.names = T, row.names = F, showNA = F, 
                 append = T)



## 4.2 Wrapper methods ####
set.seed(1234)
filterCtrl <- rfeControl(functions=rfFuncs, method="LOOCV") #, number=10)
results <- rfe(x= dfFS[,2:ncol(dfFS)],
               y= dfFS[,1], 
               sizes=ncol(dfFS)-1, 
               rfeControl=filterCtrl)
# Selecting relevant information to combine results
protsWrapper <- results$optVariables[1:50]
rfeResults <- results$variables %>%
     group_by(var) %>%
     summarize(Overall = mean(Overall)) %>%
     arrange(-Overall) %>%
  as.data.frame()
colnames(rfeResults) <- c("Proteins", "Importance")
# Saving
xlsx::write.xlsx(x = rfeResults, 
                file = "out/featureSelectionResults_LFQ.xlsx",
                sheetName = "RFE", 
                col.names = T, row.names = F, showNA = F, 
                append = T)





#...........................................................................####
# 5. Common results: top 20 ####
# Wrapper #
borutaResults <- borutaResults[1:20,]
rfeResults <- rfeResults[1:20,]
# Filter #
rocResults <- rocResults[1:20,]
gainResuls <- gainResuls[1:20,]
mrmrResults <- mrmrResults[1:20,]

## 5.1 Common elements ####
# Initial list
listVenn <- list(
  borutaResults$Proteins, 
  rfeResults$Proteins, 
  rocResults$Proteins, 
  gainResuls$Proteins,
  mrmrResults$Proteins)
names(listVenn) <- c("Boruta", "RFE", "ROC", "Information gain", "mRMR")
# Find common items
common_elements <- Reduce(intersect, listVenn)


## 5.2 Venn diagram #####
# Create graph #
paletaColor <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#1D3557") # 1
paletaColor <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E") # 2
paletaColor <- c("#08306B", "#C6612D", "#3E8E7E", "#9E9E9E", "#984EA3") # 3
paletaColor <- c("#08306B", "#C6612D", "#7FB3D5", "#D4A373", "#6C757D") # 4

vdia <- VennDiagram::venn.diagram(
  listVenn, 
  fill = paletaColor,
  filename = NULL,
  cex = 3.5,     
  cat.cex = 5,
  margin = 0.2,
  cat.dist = c(0.3, 0.3, 0.25, 0.3, 0.35)
  )

# Display graph
grid::grid.newpage()
grid::grid.draw(vdia)

# Save graph
svg("out/vennDiagram_LFQ.svg", width = 14, height = 14)
grid::grid.draw(vdia)
dev.off()


## 5.3 Tables ####
### TOP 20 proteins by method ####
listAll <- list(Boruta = borutaResults, RFE = rfeResults, ROC = rocResults, 
     Gain = gainResuls, mRMR = mrmrResults)
tableTop20 <- as.data.frame(listAll)
xlsx::write.xlsx(x = tableTop20, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Top 20", 
                 col.names = T, row.names = F, showNA = F, 
                 append = T)
###  Common proteins across all FS methods ####
listCommon <- lapply(
  listAll, 
  function(dfResult) {
    dfResult <- dfResult %>% 
      filter(Proteins %in% common_elements) %>% 
      arrange(Proteins) %>%
      dplyr::select(-Proteins)
    return(dfResult)
  }
  )
tableCommmon <- as.data.frame(listCommon)
colnames(tableCommmon) <- paste0("Importance - ", names(listAll))
colnames(tableCommmon) <- names(listAll)
xlsx::write.xlsx(x = tableCommmon, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Shared proteins", 
                 col.names = T, row.names = T, showNA = F, 
                 append = T)


## 5. 4 Heatmap #####
library(ggplot2)
# Matrix of common protein to long df
df_long <- reshape2::melt(as.matrix(tableCommmon))
colnames(df_long) <- c("Protein", "Method", "Importance")

# Scaling importance metrics to allow comparison
df_long <- df_long %>%
  group_by(Method) %>%
  mutate(Importance_scaled = (Importance - min(Importance)) / (max(Importance) - min(Importance))) %>%
  ungroup()

# Using the scaled importance to find the optimal order of the proteins based 
# on the similarity between methods (seriation package)
score_matrix_scaled <- reshape2::acast(df_long, Protein ~ Method, 
                                       value.var = "Importance_scaled")
order_rows <- seriation::get_order(seriation::seriate(score_matrix_scaled, method = "PCA"))
# Reordering 
ordered_proteins <- rownames(score_matrix_scaled)[order_rows]
df_long$Protein <- factor(df_long$Protein, levels = ordered_proteins)

# Creating heatmap with reorder rows. Values on cells are non-scaled importance. 
p1 <- ggplot2::ggplot(df_long, aes(x=Method, y=Protein, fill=Importance_scaled)) +
  ggplot2::geom_tile(color="white") +
  ggplot2::geom_text(aes(label=round(Importance, 2)), size=15) +
  ggplot2::scale_fill_gradient(
    low = "#A6CAEC", 
    high = "#08306B") +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title="",
    fill="Scaled importance") +
  xlab("") + ylab("") +
  ggplot2::theme(
    axis.text.y = element_text(face = "italic", size = 40),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 45), 
    legend.position = "bottom",
    legend.key.size = unit(4.3, "cm"),
    legend.title = element_text(size = 40),    
    legend.text = element_text(size = 40) 
    )
svg("out/heatmap_LFQ.svg", width = 12, height = 14)
print(p1)
dev.off()





#...........................................................................####
# 6. Individual predictors ####
dfFiltered <- df %>% dplyr::select(Grupos, all_of(common_elements))
cpList <- list()
for (i in colnames(dfFiltered)[-1]){
  cat("\n\n##### ", i, " {-} \n\n")
  
  dfCP <- dfFiltered[, c("Grupos", i)]
  colnames(dfCP)[2] <- "Proteina"
  cp <- cutpointr::cutpointr(data = dfCP, 
                             x  = Proteina,
                             class = Grupos, 
                             pos_class = "HF.rEF", 
                             neg_class = "HF.pEF",
                             na.rm = TRUE,
                             method = cutpointr::maximize_metric, 
                             metric = cutpointr::sum_sens_spec)
  cp$Protein <- i
  cp <- cp %>% 
    dplyr::select(Protein, acc, sensitivity, specificity, AUC, 
                  pos_class, optimal_cutpoint, direction) %>%
    rename(Accuracy = acc, 
           Sensitivity = sensitivity, 
           Specificity = specificity, 
           `Positive class` = pos_class,
           `Optimal cutpoint` = optimal_cutpoint, 
           Direction = direction) %>% 
    as.data.frame()
  
  cpList[[i]] <- cp
  cat("\n\n")
}

dfCP <- do.call(rbind, cpList)

xlsx::write.xlsx(x = dfCP, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Individual predictors", 
                 col.names = T, row.names = F, showNA = F, 
                 append = T)





#...........................................................................####
# 7. Saving the final data matrix used in the algorithms ####
xlsx::write.xlsx(x = df, 
                 file = "out/featureSelectionResults_LFQ.xlsx",
                 sheetName = "Data", 
                 col.names = T, row.names = T, showNA = F, 
                 append = T)



#...........................................................................####
# 8. Dependencies ####
pkgs <- sessionInfo()$otherPkgs
dfDependencies <- data.frame(Package = names(pkgs), Version = sapply(pkgs, function(x) x$Version))
write.table(dfDependencies, file = "dependencies.txt", append = F, sep = "\t", 
            row.names = F, quote = F, col.names = T)
















