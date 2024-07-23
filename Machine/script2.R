
work.path <- "./"; setwd(work.path) 

code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results") 
fig.path <- file.path(work.path, "Figures") 

if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

library(org.Hs.eg.db)
library(survival)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)

source(file.path(code.path, "compare.R"))

## Training Cohort ------------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -----------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

ALL_expr<- read.table(file.path(data.path, "ALL_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
ALL_surv <- read.table(file.path(data.path, "ALL_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(ALL_surv), colnames(ALL_expr))
ALL_expr <- ALL_expr[,comsam]; ALL_surv <- ALL_surv[comsam,,drop = F]

comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,])
ALL_expr <-t(ALL_expr[comgene,])

Train_set = Train_expr
Test_set =  Test_expr
ALL_set =  ALL_expr

## Public Signature ------------------------------------------------------------
pubSIG <- read.table(file.path(data.path, "public signatures.txt"), header = T)
if (!"Coef" %in% colnames(pubSIG)) pubSIG$Coef <- NA
pubSIG <- split(pubSIG[, c("SYMBOL", "Coef")], pubSIG$Model)

## My Signature ----------------------------------------------------------------
mySIGname = "MVIRGS" 
myAlgorithm = "Enet[alpha=0.6]" 

mySIG <- read.table(file.path(res.path, "fea_df.txt"), header = T) 
mySIG <- mySIG$features[mySIG$algorithm == myAlgorithm] 
mySIG <- data.frame("SYMBOL" = mySIG)

signatures <- pubSIG
signatures[[mySIGname]] <- mySIG

model <- list(); cinfo <- list() 
log.file <- file.path(res.path, "makeCox.log") 
if (file.exists(log.file)) file.remove(log.file) 
log.file <- file(log.file, open = "a")
sink(log.file, append = TRUE, type = "message")
for (i in names(signatures)){
  if (class(signatures[[i]]) == "data.frame"){
    model[[i]] <- makeCox(Features = signatures[[i]]$SYMBOL, 
                          coefs = signatures[[i]]$Coef,      
                          SIGname = i,                   
                          unmatchR = 0.2,                
                          Train_expr = Train_set,            
                          Train_surv = Train_surv,          
                          statusVar = "OS",                  
                          timeVar = "OS.time")             
  }else{
    model[[i]] = signatures[[i]]
  }
  
  cinfo[[i]] <- calCindex(model = model[[i]],                
                          name = i,                          
                          Test_expr = Test_set,             
                          Test_surv = Test_surv,            
                          Train_expr = Train_set,            
                          Train_surv = Train_surv, 
                          ALL_expr =ALL_set,          
                          ALL_surv = ALL_surv,           
                          Train_name = "TCGA_train",             
                          #Train_expr = NULL,             
                          #Train_surv = NULL,                
                          CohortVar = "Cohort",          
                          metaCohort = TRUE,              
                          statusVar = "OS",                  
                          timeVar = "OS.time")             
  message("")
}
closeAllConnections()

cinfo <- do.call(rbind, cinfo)
write.table(cinfo[,1:5], file = file.path(res.path,"cinfo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F) 
cinfo <- split(cinfo, cinfo$Cohort)

CohortCol <- c("steelblue","#CAB2D6FF",'#B8D9A9FF','brown2','blue','#FB9A99FF' ) 
names(CohortCol) <- names(cinfo)

plots <- lapply(cinfo, function(plot.data){
  plot.data$method <- 
    factor(plot.data$method,
           levels = plot.data$method[order(plot.data$C, decreasing = F)])
  
  # compares two concordance indices: the statistical test is a two-sided Student t test for dependent samples.
  C.compare <- plot.data$C[plot.data$method == mySIGname]
  se.compare <- plot.data$se[plot.data$method == mySIGname]
  n.compare <- plot.data$n[plot.data$method == mySIGname]
  RS.compare <- plot.data$RS[plot.data$method == mySIGname][[1]]
  r.combined <- unlist(lapply(plot.data$RS, function(x) cor(x, RS.compare)))
  var.combined <- plot.data$se^2 + se.compare^2 - 2*r.combined*plot.data$se*se.compare
  p <- pt(abs((plot.data$C-C.compare))/(sqrt(var.combined)), n.compare - 1, lower.tail = F) * 2
  plot.data$label <- cut(p, breaks = c(0, 0.05, 0.01, 0.001, 0.0001))
  plot.data$label <- plyr::mapvalues(x = plot.data$label,
                                     from = c("(0,0.0001]", "(0.0001,0.001]", "(0.001,0.01]", "(0.01,0.05]"), 
                                     to = c("****", "***", "**", "*"))
  
  return(ggplot(plot.data, aes(x = method, y = C, fill = Cohort)) +
           geom_errorbar(aes(ymin = C - 1.96 * se, ymax = C + 1.96 * se), width = .1) +
           geom_point(color = CohortCol[unique(plot.data$Cohort)], size = 2.5) +
           geom_text(aes(x = method, y = max(plot.data$C + 1.96 * plot.data$se - 0.05), label = label)) +
           geom_hline(yintercept = 0.6, linetype = "dashed") +
           ggtitle(label = unique(plot.data$Cohort)) +
           coord_flip() + 
           theme_classic() +
           theme(panel.border = element_rect(fill = NA, size = 1),
                 axis.title = element_blank(),
                 legend.position = "none"))
})

plot_grid(plotlist = plots, nrow = 1)
ggsave(file.path(fig.path, "comparison.pdf"), width = 12, height =14)
