
library(doParallel)
numCores=8
cl=makeCluster(numCores)
registerDoParallel(cl)

work.path <- "./"; setwd(work.path) 

code.path <- file.path(work.path, "Codes") 
data.path <- file.path(work.path, "InputData") 
res.path <- file.path(work.path, "Results") 
fig.path <- file.path(work.path, "Figures") 

if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

source(file.path(code.path, "ML.R"))

FinalModel <- c("panML", "multiCox")[2]#######

## Training Cohort ---------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因


Train_set =Train_expr
Test_set =  Test_expr

# Model training and validation -------------------------------------------
## method list --------------------------------------------------------
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 2
timeVar = "OS.time"; statusVar = "OS" 

## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed =1) 
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_expr = Train_set, 
                                 Train_surv = Train_surv, 
                                 mode = "Variable",      
                                 classVar = classVar) 
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() 
#set.seed(seed = 123) 
for (method in methods){ 
  # method <- "CoxBoost+plsRcox" 
  cat(match(method, methods), ":", method, "\n") 
  method_name = method 
  method <- strsplit(method, "\\+")[[1]] 
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2], 
                                  Train_expr = Train_expr[, selected.var], 
                                  Train_surv = Train_surv, 
                                  mode = "Model",      
                                  classVar = classVar)  
  }
  
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model, file.path(res.path, "model.rds")) 

if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) 
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) 

## Evaluate the model -----------------------------------------------------
model <- readRDS(file.path(res.path, "coxmodel.rds")) 
methodsValid <- names(model) 

RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = rbind.data.frame(Train_set,Test_set), 
                                    type = "lp") 

}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) 

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) 
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], 
                                  Test_expr = Test_set, 
                                  Test_surv = Test_surv, 
                                  Train_expr = Train_set, 
                                Train_surv = Train_surv, 
                                  Train_name = "TCGA_train", 
                                  cohortVar = "Cohort",
                                  timeVar = timeVar, 
                                  statusVar = statusVar) 
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------
Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) 
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] 
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] 

CohortCol <- c("steelblue","#CAB2D6FF",'#B8D9A9FF'
 ,'brown2','blue'
) 
names(CohortCol) <- colnames(Cindex_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, 
                    avg_Cindex = avg_Cindex, 
                    CohortCol = CohortCol, 
                    barCol = "#b15928ff", 
                    col = c("#A6CEE3FF","#FD6467"), 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 
pdf(file.path(fig.path, "heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 10, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

