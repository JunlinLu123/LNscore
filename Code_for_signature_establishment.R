library(tibble)
library(dplyr)
library(limma)
library(ggsci)
library(sva)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(glmnet)
library(MASS)
library(survival)
library(rms)
library(pROC)
library(plotROC)
library(ggthemes)
library(nomogramFormula)
library(ggrepel)
library(rstatix)
library(reshape2)
library(plyr)
library(pheatmap)
library(colorspace)
library(RImagePalette)
library(forestplot)
library(ResourceSelection)
library(ComplexHeatmap)
library(circlize)
library(rmda)
library(survminer)
library(survival)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(modEvA)
library(forestplot)

X <- intersect(
  intersect(
    intersect(
      intersect(
        intersect(
          intersect(
            intersect(rownames(mRNA_TCGA), rownames(mRNA_pair36)),
            rownames(mRNA_pair18)), 
          rownames(mRNA_GSE13507)),
        rownames(mRNA_GSE31684)), 
      rownames(mRNA_EMTAB)),
    rownames(mRNA_GSE48276)),
  rownames(mRNA_GSE52219))

X.mRNA_EMTAB <- mRNA_EMTAB[X,]
range(X.mRNA_EMTAB)
X.mRNA_GSE13507 <- mRNA_GSE13507[X,]
range(X.mRNA_GSE13507)
X.mRNA_GSE31684 <- mRNA_GSE31684[X,]
range(X.mRNA_GSE31684)
X.mRNA_GSE48276 <- mRNA_GSE48276[X,]
range(X.mRNA_GSE48276)
X.mRNA_pair42T <- mRNA_pair36[X,]
range(X.mRNA_pair42T)
X.mRNA_TCGA <- mRNA_TCGA[X,]
range(X.mRNA_TCGA)
X.mRNA_pair15 <- mRNA_pair18[X,]
X.mRNA_pair15 <- log2(2^X.mRNA_pair15+1)
range(X.mRNA_pair15)
X.mRNA_GSE52219 <- mRNA_GSE52219[X,]
range(X.mRNA_GSE52219)

combined.expr <- cbind.data.frame(X.mRNA_EMTAB,
                                  X.mRNA_GSE13507,
                                  X.mRNA_GSE31684,
                                  X.mRNA_GSE48276,
                                  X.mRNA_pair42T,
                                  X.mRNA_TCGA,
                                  X.mRNA_pair15,
                                  X.mRNA_GSE52219)


source("batchPCA.R")
mypal = pal_lancet(palette = c("lanonc"), alpha = 1)(8)
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("E_MTAB_1803","GSE13507","GSE31684","GSE48276",
                       "SYScohort1","TCGA_BLCA","SYScohort2","GSE52219"), 
                     times = c(ncol(X.mRNA_EMTAB),ncol(X.mRNA_GSE13507),ncol(X.mRNA_GSE31684),ncol(X.mRNA_GSE48276),
                               ncol(X.mRNA_pair42T),ncol(X.mRNA_TCGA),ncol(X.mRNA_pair15),ncol(X.mRNA_GSE52219))),
         fig.dir = ".",
         PCA.fig.title = "PCA before batch",
         cols = mypal,
         showID = F,
         cex = 0.7,
         showLegend = T) 


batch <- data.frame(batch = rep(c("E_MTAB_1803","GSE13507","GSE31684","GSE48276",
                                  "SYScohort1","TCGA_BLCA","SYScohort2","GSE52219"), 
                                times = c(ncol(X.mRNA_EMTAB),ncol(X.mRNA_GSE13507),ncol(X.mRNA_GSE31684),ncol(X.mRNA_GSE48276),
                                          ncol(X.mRNA_pair42T),ncol(X.mRNA_TCGA),ncol(X.mRNA_pair15),ncol(X.mRNA_GSE52219))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))


batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("E_MTAB_1803","GSE13507","GSE31684","GSE48276",
                       "SYScohort1","TCGA_BLCA","SYScohort2","GSE52219"), 
                     times = c(ncol(X.mRNA_EMTAB),ncol(X.mRNA_GSE13507),ncol(X.mRNA_GSE31684),ncol(X.mRNA_GSE48276),
                               ncol(X.mRNA_pair42T),ncol(X.mRNA_TCGA),ncol(X.mRNA_pair15),ncol(X.mRNA_GSE52219))),
         fig.dir = ".",
         PCA.fig.title = "2. PCA after batch",
         cols = mypal,
         showID = F,
         cex = 0.7,
         showLegend = T) 

combined.expr.combat<- normalizeBetweenArrays(combined.expr.combat)
write.table(combined.expr.combat,"combined.expr.combat.norm.txt",row.names = T,sep = "\t")



dim(combined.expr.combat.norm)
dim(clinical_merge)

expr.combined <- combined.expr.combat.norm_2


SampleName.pos <- clinical_merge[which(clinical_merge$lymph=="Positive"),]$patient_ID

SampleName.neg <- clinical_merge[which(clinical_merge$lymph=="Negative"),]$patient_ID

SampleName.pos.train <- sample(SampleName.pos,length(SampleName.pos)/2)
SampleName.neg.train <- sample(SampleName.neg,length(SampleName.neg)/2+0.5)
SampleName.train <- c(SampleName.pos.train,SampleName.neg.train)
expr.combined.train <- expr.combined[,SampleName.train]  #表达谱
clinical_merge.train <- clinical_merge[which(clinical_merge$patient_ID%in%SampleName.train),] #临床信息
expr.combined.train <- expr.combined.train[,clinical_merge.train$patient_ID] #调整表达谱样本顺序
clinical_merge.train$lymph <- factor(clinical_merge.train$lymph)

SampleName.pos.test <- SampleName.pos[-which(SampleName.pos%in%SampleName.pos.train)]
SampleName.neg.test <- SampleName.neg[-which(SampleName.neg%in%SampleName.neg.train)]
SampleName.test <- c(SampleName.pos.test,SampleName.neg.test)
expr.combined.test <- expr.combined[,SampleName.test]  #表达谱
clinical_merge.test <- clinical_merge[which(clinical_merge$patient_ID%in%SampleName.test),] #临床信息
expr.combined.test <- expr.combined.test[,clinical_merge.test$patient_ID] #调整表达谱样本顺序
clinical_merge.test$lymph <- factor(clinical_merge.test$lymph)



exprSet <- expr.combined.train

group_list <- as.character(clinical_merge.train$lymph)
design <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = clinical_merge.train$patient_ID

contrast.matrix <- makeContrasts("Positive-Negative",
                                 levels = design)
contrast.matrix

{
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2,coef=1,n=Inf)
  nrDEG = na.omit(tempOutput)
}

DEG <- nrDEG[which(nrDEG$P.Value<0.05&abs(nrDEG$logFC)>0.5),]


DEG.expr.train <- expr.combined.train[rownames(DEG),]
DEG.expr.test <- expr.combined.test[rownames(DEG),]

tmp.train.y <- factor(clinical_merge.train$lymph,
                      levels = c("Positive","Negative"),
                      labels= c(1,0))
tmp.train.x <- as.matrix(t(DEG.expr.train))

model.lasso <- glmnet(tmp.train.x,tmp.train.y,family="binomial",
                      nlambda = 100,alpha = 1,standardize = TRUE)
print(model.lasso)
plot(model.lasso,xvar = "lambda",label=TRUE)

cv.model <- cv.glmnet(tmp.train.x,tmp.train.y,family="binomial",
                      nlambda = 100,alpha = 1,standardize = TRUE)

plot(cv.model)
cv.model$lambda.min
coef(cv.model,s=cv.model$lambda.min)

lasso.coef <- data.frame(genename=rownames(DEG.expr.train),
                         coef=coef(cv.model,s=cv.model$lambda.min)[1:length(rownames(DEG.expr.train))])
lasso.coef <- lasso.coef[which(lasso.coef$coef!=0),]
lasso.coef
dim(lasso.coef)[1]

variable <- lasso.coef$genename

DEG.expr.test <- as.data.frame(t(DEG.expr.test))
DEG.expr.test <- rownames_to_column(DEG.expr.test,var = "patient_ID")
expr.test <- left_join(clinical_merge.test,DEG.expr.test)

clinical_merge.train
DEG.expr.variable <- DEG.expr.train[variable,]
DEG.expr.variable <- rownames_to_column(as.data.frame(t(DEG.expr.variable)),var = "patient_ID")
expr.train <- left_join(clinical_merge.train,DEG.expr.variable)
expr.train$lymph <- factor(expr.train$lymph)
expr.train$lymph <- factor(expr.train$lymph)

colnames(expr.train)

formula_1 = 
  as.formula(paste("lymph",paste(colnames(expr.train)[5:length(colnames(expr.train))],collapse="+"),sep = "~"))

model.full <- glm(formula = formula_1,
                  data=expr.train,
                  family = binomial,
                  x=T, y=T)
summary(model.full)


gene_1 <- names(which(coef(summary(model.full))[-1,4]<0.05))#获取P值有意义的基因名
formula_2 = 
  as.formula(paste("lymph",paste(gene_1,collapse="+"),sep = "~"))


model.optimal <- glm(formula = formula_2,
                     data=expr.train,
                     family = binomial,
                     x=T, y=T)

summary(model.optimal)
coef(model.optimal)


clinical_merge.train$predict <- predict(model.optimal,type="link")#计算线性评分

roc(clinical_merge.train$lymph,clinical_merge.train$predict)[c(2:3,9)]

ggplot(clinical_merge.train, aes(d = lymph, m = predict)) + 
  geom_roc(n.cuts = 0)+
  labs(x="1 - Specificity",y="sensitivity")+
  theme_few()

lrm.final <- lrm(formula = model.optimal$formula, 
                 data=expr.train,
                 x=T, y=T)
dist <- datadist(expr.train);options(datadist='dist')
nom <- nomogram(lrm.final, fun=plogis,
                fun.at=c(  .1, .25, .5, .75, .9),
                funlabel="LN+")
plot(nom,
     xfrac = .25,cex.axis = 1.5, cex.var=1.5,
     col.grid = c("red","green"))

results <- formula_lp(nomogram = nom)
clinical_merge.train$predict <- predict(model.optimal,type="link")#计算线性评分
points <- points_cal(formula = results$formula, lp = clinical_merge.train$predict)


clinical_merge.test 
expr.test 
coef(model.optimal)
clinical_merge.test$predict <- ((expr.test[,"MESP1"]*0.4106547+
                                   expr.test[,"EFEMP1"]*0.2829655-
                                   expr.test[,"PIGZ"]*0.3468013+
                                   expr.test[,"KRT23"]*0.2518427-
                                   expr.test[,"CALML3"]*0.2032473
                                 -3.0661410)*
                                  results$formula$`x^1`+results$formula$b0)

class(clinical_merge.test$predict)
roc(clinical_merge.test$lymph,clinical_merge.test$predict)[c(2:3,9)]
ggplot(clinical_merge.train, aes(d = lymph, m = predict)) + 
  geom_roc(n.cuts = 0)+
  labs(x="1 - Specificity",y="sensitivity")+
  theme_few()


up <- "#fd8d3c"
down <- "#6baed6"
  

write.csv(nrDEG,"nrDEG.csv")
colnames(nrDEG)
plot.DEG <- nrDEG[,c(1,4)]
plot.DEG <- rownames_to_column(plot.DEG,"Symbol")
plot.DEG$logPvalue <- -log10(plot.DEG$P.Value)
colnames(plot.DEG)[2] <- c("log2FC")

plot.DEG$group <- case_when(plot.DEG$log2FC > 0.5 & plot.DEG$P.Value < 0.05 ~ "Up",
                            plot.DEG$log2FC < -0.5 & plot.DEG$P.Value < 0.05 ~ "Down",
                            abs(plot.DEG$log2FC) <= 0.5 ~ "None",
                            plot.DEG$P.Value >= 0.05 ~ "None")
head(plot.DEG)

up<- filter(plot.DEG,group== "Up")
up
down<- filter(plot.DEG,group== "Down")
down
all <- rbind(up,down)

plot.DEG$size <- case_when(plot.DEG$log2FC > 0.5 & plot.DEG$P.Value < 0.05 ~ 2,
                           plot.DEG$log2FC < -0.5 & plot.DEG$P.Value < 0.05 ~ 2,
                           abs(plot.DEG$log2FC) <= 0.5 ~ 1,
                           plot.DEG$P.Value >= 0.05 ~ 1)

plot.DEG$group <- factor(plot.DEG$group,
                         levels= c( "Up", "Down", "None"),
                         ordered= T)
write.csv(plot.DEG,"差异基因.csv")

mycolor<- pal_npg()(10)
mycolor
mypal <- c(up,down,"#737373")

p1 <- ggplot(data=plot.DEG,aes(log2FC,logPvalue,color=group))+
  geom_point(size=2)+ 
  scale_colour_manual(values=mypal,
                      labels = c("Up-regulated","Down-regulated","No significance"))+
  theme_test()
p1

p2 <- p1+
  xlab(bquote(log[2](Fold~Change)))+
  ylab(bquote(-log[10](P~value)))+
  theme(legend.title=element_blank(),
        legend.position="top")
p2

p3 <- p2+
  theme(legend.text = element_text(size = 14),
        text=element_text(size = 18),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))
p3

p4 <- p3+geom_hline(yintercept = c(-log10(0.05)),
                    size= 0.8,
                    color= "gray30",
                    lty= "33",alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),
             size= 0.8,
             color= "gray30",
             lty= "33",alpha=0.5)



par(oma=c(0,1,0,0)) 
plot(model.lasso,xvar = "lambda",label=TRUE,
     cex.lab=1.3,cex.axis=1,font.lab=1,family="sans",
     col.axis="#000000",xgap.axis=0.5)

plot(cv.model,
     cex.lab=1.3,cex.axis=1,font.lab=1,family="sans",
     col.axis="#000000",xgap.axis=0.5)  


TableS1 <- as.data.frame(coef(summary(model.full)))
TableS1$OR <- exp(coef(model.full))
TableS1$LCI <- exp(confint(model.full))[,1]
TableS1$UCI <- exp(confint(model.full))[,2]
TableS1 <- TableS1[,4:7]
TableS1 <- rownames_to_column(TableS1,"Symbol")
write.csv(TableS1,"Table S1.csv")


TableS2 <- as.data.frame(coef(summary(model.optimal)))
TableS2$OR <- exp(coef(model.optimal))
TableS2$LCI <- exp(confint(model.optimal))[,1]
TableS2$UCI <- exp(confint(model.optimal))[,2]
TableS2 <- TableS2[,4:7]
TableS2 <- rownames_to_column(TableS2,"Symbol")
write.csv(TableS2,"Table S2.csv")


multi.logistic0 <- TableS1[-1,c(1,3,4,5,2)]
colnames(multi.logistic0) <- c("Gene","OR","lower.95CI","upper.95CI","p")
multi.logistic0 <- multi.logistic0[order(multi.logistic0$OR,decreasing = T),]

tabletext0 <- cbind(c("Gene",multi.logistic0$Gene),
                    c("OR",format(round(as.numeric(multi.logistic0$OR),3),nsmall = 3)),
                    c("lower 95%CI",format(round(as.numeric(multi.logistic0$lower.95CI),3),nsmall = 3)),
                    c("upper 95%CI",format(round(as.numeric(multi.logistic0$upper.95CI),3),nsmall = 3)),
                    c("pvalue",formatC(as.numeric(multi.logistic0$p), format = "e", digits = 2)))
tabletext0
forestplot(labeltext=tabletext0,
           mean=c(NA,as.numeric(multi.logistic0$OR)),
           lower=c(NA,as.numeric(multi.logistic0$lower.95CI)),
           upper=c(NA,as.numeric(multi.logistic0$upper.95CI)),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawNormalCI",
           col=fpColors(box="#41ab5d", lines="#238b45", zero = "black"),
           boxsize=0.3,#box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=1,
           lwd.zero=2,
           xticks = c(0,1,2),
           lwd.xaxis=2,
           xlab="Odds Ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "28" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(c(1,0,1.5,0), "cm")
)#7x11

multi.logistic <- TableS2[-1,c(1,3,4,5,2)]
colnames(multi.logistic) <- c("Gene","OR","lower.95CI","upper.95CI","p")

tabletext <- cbind(c("Gene",multi.logistic$Gene),
                   c("OR",format(round(as.numeric(multi.logistic$OR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(multi.logistic$lower.95CI),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(multi.logistic$upper.95CI),3),nsmall = 3)),
                   c("pvalue",formatC(as.numeric(multi.logistic$p), format = "e", digits = 2)))
tabletext
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(multi.logistic$OR)),
           lower=c(NA,as.numeric(multi.logistic$lower.95CI)),
           upper=c(NA,as.numeric(multi.logistic$upper.95CI)),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawNormalCI",
           col=fpColors(box="#41ab5d", lines="#238b45", zero = "black"),
           boxsize=0.3,#box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=1,
           lwd.zero=2,
           xticks = c(0,1,2),
           lwd.xaxis=2,
           xlab="Odds Ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "7" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1.3,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(c(1,1.5,1.5,1.5), "cm")
)#7x3


coef(model.optimal)
heatmp <- expr.train[,c("patient_ID","MESP1","EFEMP1","PIGZ","KRT23","CALML3")]
heatmp <- as.data.frame(t(column_to_rownames(heatmp,"patient_ID")))
annotation_col <- clinical_merge.train[,-4]
write.csv(annotation_col,"annotation_col.csv",row.names = F)



annotation_col <- column_to_rownames(annotation_col,"patient_ID")
colnames(annotation_col) <- c("Dataset","LN status")
annotation_col.test <- annotation_col[,2]
names(annotation_col.test) <- rownames(annotation_col)
annotation_col$Dataset <- factor(annotation_col$Dataset,
                                 levels = c("GSE13507","GSE31684","TCGA_BLCA"),
                                 labels = c("GSE13507","GSE31684","TCGA"))



heatmp<-heatmp[,rownames(annotation_col)]
range(heatmp)

write.csv(heatmp,"heatmp.csv")

dim(heatmp)
bk = unique(c(seq(-1,1, length=100)))


pal_npg("nrc")(10)
ann_colors <- list(Dataset = c(GSE13507="#984ea3",
                               GSE31684="#ff7f00",
                               TCGA="#f781bf"),
                   `LN status`=c(Negative="#4DBBD5FF",Positive="#E64B35FF"))

hcl_palettes(plot = TRUE)
sequential_hcl(5,palette = "Plasma")
qualitative_hcl(6,palette = "set2")

pheatmap(heatmp,cluster_row = FALSE, cluster_col = FALSE,
         border_color = "grey",
         annotation_col=annotation_col,
         annotation_colors = ann_colors,
         show_colnames=F,
         color = colorRampPalette(c("#1D3E64","white","#EDA200"))(100),
         fontsize = 14,
         mar=unit(c(1,1.5,1.5,0), "cm"))#6x3

colnames(clinical_score)
dim(clinical_score)
clinical_score <- clinical_score[order(clinical_score$Score,decreasing = F),]
rownames(clinical_score) <- NULL
annotation_col_2 <- clinical_score[,c(9,8,5,6,7,2,3,1,4)]
annotation_col_2 <- column_to_rownames(annotation_col_2,"patient_ID")
annotation_col_2$age <- as.numeric(annotation_col_2$age)


colnames(annotation_col_2)
col_fun0 = colorRamp2(c(0,7.5,15), c(down,"white",up))
top_annotation0 <- HeatmapAnnotation(`LN status`=annotation_col$`LN status`,
                                     dataset = annotation_col$Dataset,
                                     col=list(`LN status`=c(Positive="#a1d99b",Negative="#e5f5e0"),
                                              dataset=c(GSE13507="#efedf5",GSE31684="#dadaeb",TCGA="#bcbddc")))
ht_list =Heatmap(heatmp,
                 cluster_rows =F,
                 cluster_columns = F,
                 show_column_names = F,
                 top_annotation = top_annotation0,
                 col = col_fun0,
                 heatmap_legend_param = list(title = expression(bold(log[2]*FPKM))))
draw(ht_list,
     merge_legend = TRUE)


lrm.final <- lrm(formula = model.optimal$formula, 
                 data=expr.train,
                 x=T, y=T)
dist <- datadist(expr.train);options(datadist='dist')
nom <- nomogram(lrm.final, fun=plogis,
                fun.at=c(  .1, .25, .5, .75, .9),
                funlabel="Lymph node positive")
plot(nom,
     xfrac = .2,cex.axis = 1.2, cex.var=1.2,
     col.grid = c(up,down))


clinical_merge.train$predict <- ((expr.train[,"MESP1"]*0.4106547+
                                    expr.train[,"EFEMP1"]*0.2829655-
                                    expr.train[,"PIGZ"]*0.3468013+
                                    expr.train[,"KRT23"]*0.2518427-
                                    expr.train[,"CALML3"]*0.2032473
                                  -3.0661410)*
                                   results$formula$`x^1`+results$formula$b0)

score.train <- clinical_merge.train[order(clinical_merge.train$predict,decreasing = F),c(1,3,4)]
score.train$rank <- 1:length(rownames(score.train))
score.train$rank <- as.numeric(score.train$rank)
ggplot(score.train,aes(x=rank,y=predict,fill=lymph)) + 
  geom_bar(stat = 'identity',width = 0.7) + 
  scale_x_discrete(expand = expand_scale(mult = c(0.01,0)))+
  geom_hline(yintercept = median(score.train$predict), 
             linetype = 5, 
             size = 0.3) + 
  scale_fill_manual(values=c(down,up))+
  labs(x = "", y = "Risk Score") +
  theme_bw() +
  theme(panel.grid =element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + 
  scale_size(range=c(5,20)) +
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 14),
        text=element_text(size = 18),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))


score.test <- clinical_merge.test[order(clinical_merge.test$predict,decreasing = F),c(1,3,4)]
score.test$rank <- 1:length(rownames(score.test))
score.test$rank <- as.numeric(score.test$rank)
ggplot(score.test,aes(x=rank,y=predict,fill=lymph)) + 
  geom_bar(stat = 'identity',width = 0.7) + 
  scale_x_discrete(expand = expand_scale(mult = c(0.01,0)))+
  geom_hline(yintercept = median(score.test$predict), 
             linetype = 5, 
             size = 0.3) + 
  scale_fill_manual(values=c(down,up))+
  labs(x = "", y = "Risk Score") +
  theme_bw() +
  theme(panel.grid =element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + 
  scale_size(range=c(5,20)) +
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 14),
        text=element_text(size = 18),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm")) #8x3

fullmodel_glm <- glm(formula = model.optimal$formula, 
                     data = expr.train, 
                     family = "binomial", 
                     control = list(maxit = 50))
p.hoslem <- hoslem.test(fullmodel_glm$y, fitted(fullmodel_glm), g=10)$p.value

full_calibrate <- calibrate(lrm.final, group=expr.train$lymph) 

round(p.hoslem,3)
plot.calib.train <- data.frame(x=full_calibrate[,"predy"],
                               y=full_calibrate[,"calibrated.corrected"])
ggplot(plot.calib.train,aes(x=x,y=y))+
  geom_smooth(se=FALSE,color = up)+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8))+
  geom_abline(intercept=0,slope=1,color=down,lty=5,size=1)+
  labs(x="Nomogram predicted probability",y="Actual lymph node metastasis rate")+
  theme_few()+ 
  annotate("text",x=0.5,y=0.1,label=expression("Hosmer-Lemeshow "~italic(P)~" = "~"0.224"),size = 5) +
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 14),
        text=element_text(size = 18),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))


expr.test$lymph <- as.factor(expr.test$lymph)
fullmodel_glm_test <- glm(formula = model.optimal$formula, 
                          data = expr.test, 
                          family = "binomial", 
                          control = list(maxit = 50))
p.hoslem_test <- hoslem.test(fullmodel_glm_test$y, fitted(fullmodel_glm_test), g=10)$p.value

lrm.final.test <- lrm(formula = model.optimal$formula, 
                      data=expr.test,
                      x=T, y=T)
full_calibrate_test <- calibrate(lrm.final.test, group=expr.test$lymph) 

round(p.hoslem_test,3)
plot.calib.test <- data.frame(x=full_calibrate_test[,"predy"],
                              y=full_calibrate_test[,"calibrated.corrected"])
ggplot(plot.calib.test,aes(x=x,y=y))+
  geom_smooth(se=FALSE,color = up)+
  scale_x_continuous(breaks = c(0.2,0.4,0.6))+
  geom_abline(intercept=0,slope=1,color=down,lty=5,size=1)+
  labs(x="Nomogram predicted probability",y="Actual lymph node metastasis rate")+
  theme_few()+ 
  annotate("text",x=0.4,y=0.1,label=expression("Hosmer-Lemeshow "~italic(P)~" = "~"0.178"),size = 5) +
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 14),
        text=element_text(size = 18),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))#6x5.5

roc(clinical_merge.train$lymph,clinical_merge.train$predict)[c(2:3,9)]
ggplot(clinical_merge.train, aes(d = lymph, m = predict)) + 
  geom_roc(n.cuts = 0,size=0.8,color=up)+
  geom_abline(intercept=0,slope=1 )+
  labs(x="1 - Specificity",y="Sensitivity")+
  theme_few()+
  annotate("text",x=0.5,y=0.1,
           label="AUC = 0.7811",
           size=4.5)+
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 12),
        text=element_text(size = 16),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))#5x4.5

roc(clinical_merge.test$lymph,clinical_merge.test$predict)[c(2:3,9)]
ggplot(clinical_merge.test, aes(d = lymph, m = predict)) + 
  geom_roc(n.cuts = 0,size=0.8,color=up)+
  geom_abline(intercept=0,slope=1 )+
  labs(x="1 - Specificity",y="Sensitivity")+
  theme_few()+
  annotate("text",x=0.5,y=0.1,
           label="AUC = 0.6611",
           size=4.5)+
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 12),
        text=element_text(size = 16),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))

data <- clinical_merge.train
data$lymph[which(data$lymph=="Positive")]=1
data$lymph[which(data$lymph=="Negative")]=0
data$lymph <- as.numeric(data$lymph)
Model.train <- decision_curve(lymph ~ predict,
                              data=data,
                              family = binomial(link = 'logit'),
                              thresholds = seq(0,1,by = 0.01),
                              confidence.intervals = 0.95,
                              study.design = 'case-control',
                              population.prevalence = 0.25)

plot_decision_curve(list(Model),
                    curve.names = c('Nomogram'),
                    cost.benefit.axis = F,
                    confidence.intervals = F,
                    standardize = F,
                    col = c("#9e9ac8"),
                    lty = 1, 
                    xlab = 'Risk Threshold',
                    lwd	= 2,
                    font.axis = 1,
                    legend.position = "topright")



clinical_merge.train$type <- "train"
clinical_merge.test$type <- "test"
clinical_score <- rbind(clinical_merge.train,clinical_merge.test)
write.csv(clinical_score,"clinical_score.csv",row.names = F)

colnames(Sec_Figure6)
Sec_Figure6$Image <- as.factor(Sec_Figure6$Image)
sec_fig6_LNp <- Sec_Figure6[which(Sec_Figure6$LN_c=="positive"),]
sec_fig6_LNp$rank <- c(1:length(sec_fig6_LNp$rank))
sec_fig6_LNn <- Sec_Figure6[which(Sec_Figure6$LN_c=="negative"),]
sec_fig6_LNn$rank <- c(1:length(sec_fig6_LNn$rank))
sec_fig6_LN <- rbind(sec_fig6_LNp,sec_fig6_LNn)

sec_fig6_LN$score <- sec_fig6_LN$score+5.1
ggplot(sec_fig6_LN,aes(x=rank,y=score,fill=Image)) + 
  geom_bar(stat = 'identity',width = 0.9)+
  scale_fill_manual(values = c("#41b6c4","#fc8d59"))+
  scale_x_discrete(expand = expansion(mult = c(0.01,0)))+
  geom_hline(yintercept = c(10.5), 
             linetype = 5, 
             size = 0.5,
             alpha = 0.5)+
  labs(fill="CT/MRI\nDiagnosis",
       x = "", y = "Lymph node score")+
  theme_bw() +
  theme(panel.grid =element_blank())+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = 
                                           rgb(1,1,1,alpha=0.001)),
        axis.title =element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))+
  facet_grid(~LN_c,scales = "free_x",space="free_x",shrink=F)+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0.1,5.1,10.1,15.1),
                     labels = c("-5","0","5","10"))+
  #theme(legend.position=c(0.02,0.65),legend.justification = c(0,0),
  #     legend.box = "horizontal")+
  guides(fill = guide_legend(reverse=TRUE))#3x12


colnames(Sec_Figure6)
Sec_Figure6$Image <- as.factor(Sec_Figure6$Image)
sec_fig6_LNp <- Sec_Figure6[which(Sec_Figure6$LN_c=="positive"),]
sec_fig6_LNp$rank <- c(1:length(sec_fig6_LNp$rank))
sec_fig6_LNn <- Sec_Figure6[which(Sec_Figure6$LN_c=="negative"),]
sec_fig6_LNn$rank <- c(1:length(sec_fig6_LNn$rank))
sec_fig6_LN <- rbind(sec_fig6_LNp,sec_fig6_LNn)

ggplot(sec_fig6_LNp,aes(x=rank,y=score,fill=Image)) + 
  geom_bar(stat = 'identity',width = 0.9)+
  scale_fill_manual(values = c("#41b6c4","#fc8d59"))+
  scale_x_discrete(expand = expansion(mult = c(0.1,0.1)))+
  geom_hline(yintercept = 5.5, 
             linetype = 5, 
             size = 0.5,
             alpha = 0.5)+
  labs(fill="CT/MRI diagnosis",
       x = "", y = "Lymph node gene score")+
  theme_bw() +
  theme(panel.grid =element_blank())+
  theme(legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = 
                                           rgb(1,1,1,alpha=0.001)),
        axis.title =element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,5,10))+
  guides(fill = guide_legend(reverse=TRUE))#4x3

colnames(Sec_Figure6)
Sec_Figure6$Image <- as.factor(Sec_Figure6$Image)
sec_fig6_LNp <- Sec_Figure6[which(Sec_Figure6$LN_c=="positive"),]
sec_fig6_LNp$rank <- c(1:length(sec_fig6_LNp$rank))
sec_fig6_LNn <- Sec_Figure6[which(Sec_Figure6$LN_c=="negative"),]
sec_fig6_LNn$rank <- c(1:length(sec_fig6_LNn$rank))

sec_fig6_LNn$score <- sec_fig6_LNn$score+5.1

ggplot(sec_fig6_LNn,aes(x=rank,y=score,fill=Image)) + 
  geom_bar(stat = 'identity',width = 0.9)+
  scale_fill_manual(values = c("#41b6c4","#fc8d59"))+
  geom_hline(yintercept = 10.6, 
             linetype = 5, 
             size = 0.4,
             alpha = 0.5)+
  scale_x_discrete(expand = expansion(mult = c(0.01,0.01)))+
  labs(fill="CT/MRI diagnosis",
       x = "", y = "Lymph node gene score")+
  theme_bw() +
  theme(panel.grid =element_blank())+
  theme(legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = 
                                           rgb(1,1,1,alpha=0.001)),
        axis.title =element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())+
  scale_y_continuous(limits = c(0,15),
                     breaks = c(0.1,5.1,10.1,15.1),
                     labels = c("-5","0","5","10"))+
  guides(fill = guide_legend(reverse=TRUE))#4x10

colnames(clinical_score)
clinical_score_LNn <- clinical_score[which(clinical_score$lymph=="Negative"),]
clinical_score_total_LNn <- clinical_score_total[
  which(clinical_score_total$patient_ID%in%clinical_score_LNn$patient_ID),]
rownames(clinical_score_total_LNn) <- NULL

ggsurvplot(survfit(Surv(OS.time,OS)~Score_c,
                   clinical_score_total_LNn),
           pval = T,
           pval.size = 5,
           conf.int = F,
           ylab = "Overall survival",
           xlab = "Time in days",
           ggtheme = theme_bw(),
           legend.title = NULL,
           legend.labs = c("High score","Low score"),
           font.legend = c(12),
           font.x = c(14),
           font.y = c(14),
           font.tickslab = c(12),
           risk.table = T)#6x5

colnames(clinical_score_total_LNn)
rownames(clinical_score_total_LNn) <- NULL
clinical_score_total_LNn_TCGA <- clinical_score_total_LNn[c(1:236),]
ggsurvplot(survfit(Surv(OS.time,OS)~Score_c,
                   clinical_score_total_LNn_TCGA),
           pval = T,
           pval.size = 5,
           conf.int = F,
           ylab = "Overall survival",
           xlab = "Time in days",
           ggtheme = theme_bw(),
           legend.title = NULL,
           legend.labs = c("High score","Low score"),
           font.legend = c(12),
           font.x = c(14),
           font.y = c(14),
           font.tickslab = c(12),
           risk.table = T)#6x5

clinical_score_total_LNn_GSE13507 <- clinical_score_total_LNn[c(237:385),]
ggsurvplot(survfit(Surv(OS.time,OS)~Score_c,
                   clinical_score_total_LNn_GSE13507),
           pval = T,
           pval.size = 5,
           conf.int = F,
           ylab = "Overall survival",
           xlab = "Time in days",
           ggtheme = theme_bw(),
           legend.title = NULL,
           legend.labs = c("High score","Low score"),
           font.legend = c(12),
           font.x = c(14),
           font.y = c(14),
           font.tickslab = c(12),
           risk.table = T)#6x5

clinical_score_total_LNn_GSE31684 <- clinical_score_total_LNn[c(386:434),]
ggsurvplot(survfit(Surv(OS.time,OS)~Score_c,
                   clinical_score_total_LNn_GSE31684),
           pval = T,
           pval.size = 5,
           conf.int = F,
           ylab = "Overall survival",
           xlab = "Time in days",
           ggtheme = theme_bw(),
           legend.title = NULL,
           legend.labs = c("High score","Low score"),
           font.legend = c(12),
           font.x = c(14),
           font.y = c(14),
           font.tickslab = c(12),
           risk.table = T)#6x5

colnames(F6survival)
ggsurvplot(survfit(Surv(OS.time,OS)~score_c,
                   F6survival),
           pval = T,
           pval.size = 5,
           conf.int = F,
           ylab = "Overall survival",
           xlab = "Time in days",
           ggtheme = theme_bw(),
           legend.title = NULL,
           legend.labs = c("High score","Low score"),
           font.legend = c(12),
           font.x = c(14),
           font.y = c(14),
           font.tickslab = c(12),
           risk.table = T)#6x5


#绘制PR曲线-SYSMH
sec_pr_SYS <- sec_fig6_LN
sec_pr_SYS$LN_c[which(sec_pr_SYS$LN_c=="positive")] <- 1
sec_pr_SYS$LN_c[which(sec_pr_SYS$LN_c=="negative")] <- 0

sec_pr_SYS$LN_c <- as.numeric(sec_pr_SYS$LN_c)
class(sec_pr_SYS$LN_c)
model.PR <- glm(formula = LN_c~score,
                data=sec_pr_SYS,
                family = binomial,
                x=T, y=T)
pred <- predict(model.PR,type= "response")
pr_AUC <- AUC(obs = sec_pr_SYS$LN_c,
              pred = as.numeric(pred),
              curve = "PR",
              main="PR curve")
pr_AUC_data <- data.frame(precision = pr_AUC$thresholds$precision,
                          recall = pr_AUC$thresholds$sensitivity,
                          prauc = pr_AUC$AUC)
mod_pr_AUC_data <- pr_AUC_data[c(nrow(pr_AUC_data):1),]
#画图
ggplot(mod_pr_AUC_data,aes(recall,precision))+
  geom_line(color="#fc8d59",
            size=1)+
  theme_few()+
  ylim(0,1)+
  xlim(0,1)+
  geom_abline(intercept=1,slope=-1 )+
  labs(y="Precision",x="Recall")+
  annotate("text",x=0.5,y=0.1,
           label="PR-AUC = 0.774",
           size=4.5)+
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 12),
        text=element_text(size = 16),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))#5x4.5

#绘制PR曲线-Training
sec_pr_train <- clinical_merge.train
sec_pr_train$lymph[which(sec_pr_train$lymph=="Positive")] <- 1
sec_pr_train$lymph[which(sec_pr_train$lymph=="Negative")] <- 0
sec_pr_train$lymph <- as.numeric(sec_pr_train$lymph)
model.PR <- glm(formula = lymph~predict,
                data=sec_pr_train,
                family = binomial,
                x=T, y=T)
pred <- predict(model.PR,type= "response")
pr_AUC <- AUC(obs = sec_pr_train$lymph,
              pred = as.numeric(pred),
              curve = "PR",
              main="PR curve")
pr_AUC_data <- data.frame(precision = pr_AUC$thresholds$precision,
                          recall = pr_AUC$thresholds$sensitivity,
                          prauc = pr_AUC$AUC)
mod_pr_AUC_data <- pr_AUC_data[c(nrow(pr_AUC_data):1),]
ggplot(mod_pr_AUC_data,aes(recall,precision))+
  geom_line(color="#fc8d59",
            size=1)+
  theme_few()+
  ylim(0,1)+
  xlim(0,1)+
  geom_abline(intercept=1,slope=-1 )+
  labs(y="Precision",x="Recall")+
  annotate("text",x=0.5,y=0.1,
           label="PR-AUC = 0.774",
           size=4.5)+
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 12),
        text=element_text(size = 16),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))#5x4.5

sec_pr_test <- clinical_merge.test
sec_pr_test$lymph[which(sec_pr_test$lymph=="Positive")] <- 1
sec_pr_test$lymph[which(sec_pr_test$lymph=="Negative")] <- 0
sec_pr_test$lymph <- as.numeric(sec_pr_test$lymph)
model.PR <- glm(formula = lymph~predict,
                data=sec_pr_test,
                family = binomial,
                x=T, y=T)
pred <- predict(model.PR,type= "response")
pr_AUC <- AUC(obs = sec_pr_test$lymph,
              pred = as.numeric(pred),
              curve = "PR",
              main="PR curve")
pr_AUC_data <- data.frame(precision = pr_AUC$thresholds$precision,
                          recall = pr_AUC$thresholds$sensitivity,
                          prauc = pr_AUC$AUC)
mod_pr_AUC_data <- pr_AUC_data[c(nrow(pr_AUC_data):1),]
ggplot(mod_pr_AUC_data,aes(recall,precision))+
  geom_line(color="#fc8d59",
            size=1)+
  theme_few()+
  ylim(0,1)+
  xlim(0,1)+
  geom_abline(intercept=1,slope=-1 )+
  labs(y="Precision",x="Recall")+
  annotate("text",x=0.5,y=0.1,
           label="PR-AUC = 0.774",
           size=4.5)+
  theme(legend.position=c(0.02,0.8),legend.justification = c(0,0),
        legend.text = element_text(size = 12),
        text=element_text(size = 16),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(1.5,units = "mm"))#5x4.5


dim(combined.expr.combat.norm_2)
signature.expr <- rownames_to_column(combined.expr.combat.norm_2,"Gene")
signature.expr <- left_join(othersignature,signature.expr,"Gene")
dim(signature.expr)
signature.expr <- column_to_rownames(signature.expr,"Gene")
signature.expr <- as.data.frame(t(signature.expr))
signature.expr <- rownames_to_column(signature.expr,"patient_ID")
signature_clinical <- left_join(clinical_score[,c(1,3)],signature.expr)
signature_clinical$lymph[which(signature_clinical$lymph=="Positive")] <- 1
signature_clinical$lymph[which(signature_clinical$lymph=="Negative")] <- 0
signature_clinical$lymph <- as.numeric(signature_clinical$lymph)

gene_collection <- indi_gene$A45      
gene_collection <- gene_collection[complete.cases(gene_collection)]

formula_sig = 
  as.formula(paste("lymph",
                   paste(gene_collection,collapse="+"),
                   sep = "~"))
formula_sig
model.signature <- glm(formula = formula_sig,
                       data=signature_clinical,
                       family = binomial,
                       x=T, y=T)
roc(signature_clinical$lymph,
    predict(model.signature,type="link"))[c(9)]
roc <- roc(signature_clinical$lymph,
           predict(model.signature,type="link"))
ci(roc)

signature_ROC
multi.logistic0 <- TableS1[-1,c(1,3,4,5,2)]
colnames(multi.logistic0) <- c("Gene","OR","lower.95CI","upper.95CI","p")
multi.logistic0 <- multi.logistic0[order(multi.logistic0$OR,decreasing = T),]
signature_ROC_1<-signature_ROC[1:27,]
signature_ROC_2<-signature_ROC[28:55,]
dim(signature_ROC)

tabletext_ROC_1 <- cbind(c("Signature",signature_ROC_1$Name),
                         c("GN",signature_ROC_1$No.),
                         c("AUC",format(round(as.numeric(signature_ROC_1$ROC),3),nsmall = 3)),
                         c("L95%CI",format(round(as.numeric(signature_ROC_1$lower),3),nsmall = 3)),
                         c("U95%CI",format(round(as.numeric(signature_ROC_1$higher),3),nsmall = 3)))
tabletext_ROC_1
forestplot(labeltext=tabletext_ROC_1,
           mean=c(NA,as.numeric(signature_ROC_1$ROC)),
           lower=c(NA,as.numeric(signature_ROC_1$lower)),
           upper=c(NA,as.numeric(signature_ROC_1$higher)),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawNormalCI",
           col=fpColors(box="#41ab5d", lines="#238b45", zero = "black"),
           boxsize=0.3,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0.7,
           lwd.zero=2,
           xticks = c(0.5,0.8),
           lwd.xaxis=2,
           xlab="Odds Ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "29" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(c(1,0,1.5,0), "cm")
)#12x9

tabletext_ROC_2 <- cbind(c("Signature",signature_ROC_2$Name),
                         c("GN",signature_ROC_2$No.),
                         c("AUC",format(round(as.numeric(signature_ROC_2$ROC),3),nsmall = 3)),
                         c("L95%CI",format(round(as.numeric(signature_ROC_2$lower),3),nsmall = 3)),
                         c("U95%CI",format(round(as.numeric(signature_ROC_2$higher),3),nsmall = 3)))
tabletext_ROC_2
forestplot(labeltext=tabletext_ROC_2,
           mean=c(NA,as.numeric(signature_ROC_2$ROC)),
           lower=c(NA,as.numeric(signature_ROC_2$lower)),
           upper=c(NA,as.numeric(signature_ROC_2$higher)),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawNormalCI",
           col=fpColors(box="#41ab5d", lines="#238b45", zero = "black"),
           boxsize=0.3,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0.7,
           lwd.zero=2,
           xticks = c(0.5,0.8),
           lwd.xaxis=2,
           xlab="Odds Ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "30" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(1,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(c(1,0,1.5,0), "cm")
)#13x9








