if(FALSE){
  BiocManager::install(version = '3.13')
  #BiocManager::install("reactome.db")
  BiocManager::install("DESeq2")

  library(devtools)
  install_github("ctlab/fgsea")
  
  BiocManager::install("ReactomeGSA")
  BiocManager::install('EnhancedVolcano')
}

library(DESeq2)
library(sqldf)
library(stringr)
library(umap)
library(gplots)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)

map_genotype_human <- list(n2="WT",g="shk-1 KO")


# Wild type, no infection
# op_n2_5*
# Mutant, no infection
# op_803_1*
# Wild type, infection
# vb_n2_1 *
# Mutant, infection
# vb_803_1*



###########################################################################################
################################## Read data ##############################################
###########################################################################################


cnt <- read.table("count.csv",sep="\t", stringsAsFactors = FALSE)
rownames(cnt) <- cnt[,1]
colnames(cnt) <- str_replace_all(cnt[1,],"_1.fq.gzAligned.sortedByCoord.out.bam","")
cnt <- cnt[-1,-c(1:6)]
write.table(cnt ,"cnt_clean.csv")
cnt <- read.table("cnt_clean.csv", stringsAsFactors = FALSE)

cond <- read.csv("cond.csv",stringsAsFactors = FALSE)[,-1]
rownames(cond) <- colnames(cnt)
#cond$samp
#rownames(cond) == colnames(cnt)

colSums(cnt)/1e6
# all samples have good depth

## Kick out CC81-3
# keep <- cond$sample!="CC81_3"
# cond <- cond[keep,]
# cnt <- cnt[,keep]

write.table(cond ,"cond_clean.csv")

cond$totalcond <- factor(sprintf("%s_%s",cond$strain, cond$infected))

###########################################################################################
################################## ID mapping #############################################
###########################################################################################


gid <- read.csv("c_elegans.PRJNA13758.current.geneIDs.txt", header = FALSE)
gid <- gid[gid[,5]=="Live" ,c(2,3,4,6)]  
colnames(gid) <- c("Gene.stable.ID","genesym","transcript","genetype")
gid_entrez <- read.csv("mart_export.txt",sep="\t")[,c(2,3)]
gid <- unique(merge(gid_entrez,gid))
gid

# use genesym, otherwise gene stable id, otherwise wormbaseid, for counts
gid <- gid[gid$Gene.stable.ID %in% rownames(cnt),]
gid <- gid[!duplicated(gid$transcript),]
rownames(gid) <- gid$Gene.stable.ID
gid$genesym[gid$genesym==""] <- gid$transcript[gid$genesym==""]

gid <- gid[!duplicated(gid$genesym),]

newnames <- gid[rownames(cnt),]$genesym
newnames[is.na(newnames)] <- rownames(cnt)[is.na(newnames)]
df_name_mapping <- data.frame(wbid=rownames(cnt), genesym=newnames)
rownames(cnt) <- newnames

protein_coding <- c(gid$Gene.stable.ID[gid$genetype=="protein_coding_gene"],
                    gid$genesym[gid$genetype=="protein_coding_gene"])




###########################################################################################
################################## pairwise DE ############################################
###########################################################################################


dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cond, design = ~totalcond)
dds <- DESeq(dds)

treatments <- unique(cond$totalcond)
for(i in 1:length(treatments)){
  for(j in 1:length(treatments)){
    if(i!=j){
      t_i <- as.character(treatments[i])
      t_j <- as.character(treatments[j])
      print(c(t_i, t_j))
      r <- as.data.frame(results(dds, contrast = c("totalcond",t_i,t_j)))  #no longer ... #reversed as requested ... should I?
      r <- r[order(r$pvalue),]

      r <- r[rownames(r) %in% protein_coding,]

      # volcano plot      
      pdf(sprintf("out.de/%s___vs___%s.pdf",t_i,t_j))
      plot(r$log2FoldChange, -log10(r$padj), pch=19,cex=0.5,
           xlab="log2(FC)",ylab="-log10(p.adj)",
           main=sprintf("%s vs %s",t_i,t_j))
      dev.off()
      
      r <- r[!is.na(r$pvalue),]
      r$baseMean <- format(r$baseMean, scientific=T, digits=3)  # sprintf(r$baseMean, fmt = '%#.2f')
      r$log2FoldChange <-  sprintf(r$log2FoldChange, fmt = '%#.2f')
      r$padj <- format(r$padj, scientific=T, digits=3)
      r$pvalue <- format(r$pvalue, scientific=T, digits=3)
      write.csv(r[,c("baseMean","log2FoldChange","pvalue","padj")], sprintf("out.de/%s___vs___%s.csv",t_i,t_j))
    }
  }
}



################################################################################
# volcano: inf/uninf
################################################################################

pcutoff <- 1e-20

set.seed(666)
allplots <- list()
for(curg in c("g","n2")){
  res <- read.csv(sprintf("out.de/%s_inf___vs___%s_uninf.csv", curg, curg))
  
  res$highlight <- ""
  res$highlight[res$X %in% c("dod-19", "dod-22", "dod-24", "clec-85", "clec-4", "clec-66", "clec-67", "irg-4", "irg-5", "fil-1", "cpr-3", "lys-2", "lys-1", 
                             "fat-6", "fat-7", "acdh-1", "acdh-2", "ech-6", "asm-3")] <- "target"
  
  print(res[res$X %in% "irg-4",])

  res$color <- "_"
  res$color[res$pvalue<pcutoff & res$log2FoldChange<0] <- "down"
  res$color[res$pvalue<pcutoff & res$log2FoldChange>0] <- "up"
  res$color[res$highlight!=""] <- "xhighlight"
  res$color <- factor(res$color, levels = c("_","up","down","xhighlight"))
  res <- res[order(res$color),]
  
  res$logp <- -log10(res$pvalue)
  maxrange <- max(res$logp[!is.infinite(res$logp)], na.rm=TRUE)
  print(maxrange)
  res$logp[is.infinite(res$logp)] <- maxrange + runif(sum(is.infinite(res$logp)),0,50)

  onep <- ggplot(res, aes(log2FoldChange, logp, color=color, label=X)) + geom_point() + 
    theme_bw() + xlim(c(-7,7)) + ylim(0,370)+
    xlab(sprintf("Log2 FC, %s inf/uninf",map_genotype_human[[curg]])) +
    scale_color_manual(values = c("_"="#DDDDDD","down"="#CCCC00","up"="lightblue","xhighlight"="black"))+
    ggrepel::geom_text_repel(data=res[res$highlight!="",],color="black",min.segment.length = unit(0, 'lines'), max.overlaps = 100) +  
    theme(legend.position = "none")
  allplots[[paste("1",curg)]] <- onep

  print(paste("========",sprintf("out.de/%s_inf___vs___%s_uninf.csv", curg, curg)))
  print(paste("all DE#",sum(res$pvalue < pcutoff)))
  print(paste("up in inf",sum(res$pvalue < pcutoff & res$log2FoldChange>0)))    
  print(paste("down in inf",sum(res$pvalue < pcutoff & res$log2FoldChange<0)))
  
}

for(curcond in c("inf","uninf")){
  res <- read.csv(sprintf("out.de/g_%s___vs___n2_%s.csv", curcond, curcond))
  
  res$highlight <- ""
  
  #res$highlight[res$X %in% c("zip-8","cav-1",sprintf("catp-%s",1:8))] <- "target"
  res$highlight[res$X %in% c("zip-8","cav-1","fgt-1")] <- "target"
  
  
  res$color <- "_"
  res$color[res$pvalue<pcutoff & res$log2FoldChange<0] <- "down"
  res$color[res$pvalue<pcutoff & res$log2FoldChange>0] <- "up"
  res$color[res$highlight!=""] <- "xhighlight"
  res$color <- factor(res$color, levels = c("_","up","down","xhighlight"))
  res <- res[order(res$color),]
  
  onep <- ggplot(res, aes(log2FoldChange, -log10(pvalue), color=color, label=X)) + geom_point() + 
    theme_bw() + xlim(c(-7,7)) + #ylim(0,300)+
    xlab(sprintf("Log2 FC, shk-1 KO %s/WT %s",curcond,curcond)) +
    scale_color_manual(values = c("_"="#DDDDDD","down"="#CCCC00","up"="lightblue","xhighlight"="black"))+
    ggrepel::geom_text_repel(data=res[res$highlight!="",],color="black",min.segment.length = unit(0, 'lines'), max.overlaps = 100) +  #
    theme(legend.position = "none")
  allplots[[paste("2",curcond)]] <- onep

  print(paste("========",sprintf("out.de/g_%s___vs___n2_%s.csv", curcond, curcond)))
  print(paste("all DE#",sum(res$pvalue < pcutoff)))
  print(paste("up in KO",sum(res$pvalue < pcutoff & res$log2FoldChange>0)))    
  print(paste("down in KO",sum(res$pvalue < pcutoff & res$log2FoldChange<0)))
  
}
ptot <- egg::ggarrange(plots=allplots, ncol=2)
ggsave("newout/volcano_pickedgenes.pdf",ptot, width = 8,height = 8)




################################################################################
#Venn diagram: infection induced DE genes in WT vs shk-1.
################################################################################


res_ko <- read.csv(sprintf("out.de/g_inf___vs___g_uninf.csv"))
res_wt <- read.csv(sprintf("out.de/n2_inf___vs___n2_uninf.csv"))

list4venn_down <- list(
  ko=res_ko$X[res_ko$pvalue < pcutoff & res_ko$log2FoldChange<0],
  wt_inf_uninf=res_wt$X[res_wt$pvalue < pcutoff & res_wt$log2FoldChange<0]
  )
list4venn_up <- list(
  ko=res_ko$X[res_ko$pvalue < pcutoff & res_ko$log2FoldChange>0],
  wt_inf_uninf=res_wt$X[res_wt$pvalue < pcutoff & res_wt$log2FoldChange>0]
)
ptot <- (ggVennDiagram(list4venn_down) + ggtitle("Down")) | (ggVennDiagram(list4venn_up) + ggtitle("Up"))
ptot 
ggsave(plot = ptot, "newout/venn_inf_vs_uninf.pdf")



res_uninf <- read.csv(sprintf("out.de/g_uninf___vs___n2_uninf.csv"))
res_inf <- read.csv(sprintf("out.de/g_inf___vs___n2_inf.csv"))

list4venn_down <- list(
  uninf=res_uninf$X[res_uninf$pvalue < pcutoff & res_uninf$log2FoldChange<0],
  inf=res_inf$X[res_inf$pvalue < pcutoff & res_inf$log2FoldChange<0]
)
list4venn_up <- list(
  uninf=res_uninf$X[res_uninf$pvalue < pcutoff & res_uninf$log2FoldChange>0],
  inf=res_inf$X[res_inf$pvalue < pcutoff & res_inf$log2FoldChange>0]
)
ptot <- (ggVennDiagram(list4venn_down) + ggtitle("Down")) | (ggVennDiagram(list4venn_up) + ggtitle("Up"))
ptot 
ggsave(plot = ptot, "newout/venn_ko_vs_wt.pdf")






################################################################################
# Comparison of differentially expressed genes by infection in wt and mutants
################################################################################

keep <- cond$strain=="g"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~infected)
dds <- DESeq(dds)
coef(dds)
resLFC.g <- lfcShrink(dds, coef = c("infected_uninf_vs_inf"), type="apeglm") #cite https://doi.org/10.1093/bioinformatics/bty895

keep <- cond$strain=="n2"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~infected)
dds <- DESeq(dds)
coef(dds)
resLFC.wt <- lfcShrink(dds, coef = c("infected_uninf_vs_inf"), type="apeglm") #cite https://doi.org/10.1093/bioinformatics/bty895

df <- data.frame(
  gene=rownames(resLFC.g),
  ko=-resLFC.g$log2FoldChange,  #upregulated now positive
  wt=-resLFC.wt$log2FoldChange,
  p_ko = resLFC.g$pvalue,
  p_wt = resLFC.wt$pvalue
)
#ggplot(df[df$p_ko<pcutoff | df$p_wt<pcutoff,], aes(ko,wt)) + geom_point()


df$class <- "_"
#df$class[df$gene %in% c("fat-6", "fat-7", "acdh-1", "acdh-2", "ech-6", "mtl-2")] <- "highlight"
#df$class[str_starts(df$gene,"dod-")] <- "dod"
#df$class[str_starts(df$gene,"irg-")] <- "irg"
#df$class[str_starts(df$gene,"clec-")] <- "clec"
#df$class[str_starts(df$gene,"lys-")] <- "lys"


df$class[df$gene %in% c("clec-85", "clec-4", "clec-66", "dod-19", "dod-24", "dod-22", 
                        "irg-6", "irg-4", "lys-3", "lys-2", "fat-6", "fat-7", "acdh-1", "acdh-2", "ech-6", "mtl-2")] <- "highlight"

df <- df[order(df$class),]
ggplot(df[df$p_ko<pcutoff | df$p_wt<pcutoff,], aes(ko,wt)) + 
  geom_point(aes(ko,wt,color=class, label=gene)) + 
  geom_smooth(df[df$class=="highlight",],method = "lm", mapping = aes(color=class), color="blue")+
  geom_abline(slope=1)+
  ggrepel::geom_text_repel(mapping = aes(ko,wt,color=class, label=gene), data=df[df$class=="highlight",], color="red") +
  scale_color_manual(values = c(
    "_" = "#DDDDDD",
    "highlight" = "black",
    "dod" = "blue",
    "irg" = "green",
    "clec" = "gold",
    "lys" = "darkcyan"
  )) + 
  xlim(c(-5,12))+
  ylim(c(-5,12))+
  xlab("shk-1 KO log2(FC)") + 
  ylab("WT log2(FC)") +
  theme_bw() +
  theme(legend.position="none")
ggsave("newout/fig5.pdf", width = 4, height = 4)



################################################################################
#   Volcano plot, highlighting zip-8, cav-1, catp-1, catp-3, catp-4, and catp-7 w/ or w/o infection-
################################################################################

if(FALSE){
  glist <- c("zip-8", "cav-1","fgt-1")

  tab <- read.csv("out.de/g_inf___vs___g_uninf.csv")
  tab$pvalue[tab$pvalue==0 | tab$pvalue<1e-300] <- 1e-300
  tab$col <- "black"
  tab$col[tab$X %in% glist] <- "blue"
  tab <- tab[order(tab$col),]
  p1 <- ggplot(tab, aes(log2FoldChange, -log10(pvalue), color=col, label=X)) + 
    scale_color_manual(values=c("gray", "black")) +
    geom_point() + 
    ggrepel::geom_text_repel(data = tab[tab$X %in% glist,]) + 
    xlab("log2 FC, shk-1 KO") + ylab("-log10 (p.val), shk-1 KO")+
    xlim(-6,13)
  
  
  tab <- read.csv("out.de/n2_inf___vs___n2_uninf.csv")
  tab$pvalue[tab$pvalue==0 | tab$pvalue<1e-300] <- 1e-300
  tab$col <- "black"
  tab$col[tab$X %in% glist] <- "blue"
  tab <- tab[order(tab$col),]
  p2 <- ggplot(tab, aes(log2FoldChange, -log10(pvalue), color=col, label=X)) + 
    scale_color_manual(values=c("gray", "black")) +
    geom_point() + 
    ggrepel::geom_text_repel(data = tab[tab$X %in% glist,]) + 
    xlab("log2 FC, WT") + ylab("-log10 (p.val), WT")+
    xlim(-6,13)
  
  ptot <- egg::ggarrange(p1,p2)
  ggsave("out.plot/volcano_inf.pdf", ptot, width = 8, height = 10)  
}





###########################################################################################
############################### GO analysis for each strain ############################### 
###########################################################################################



library(enrichR)


setEnrichrSite("WormEnrichr")
listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","GO_Molecular_Function_2018") #"GO_Biological_Process_2018",#,"KEGG_2019",

list_comp <-  c(
  "g_uninf___vs___n2_uninf.csv",
  "g_inf___vs___n2_inf.csv",
  
  "g_inf___vs___g_uninf.csv",
  "n2_inf___vs___n2_uninf.csv"
)

newoutgo <- NULL
allGOdf <- NULL
for(one_comp in list_comp){
  
  
  print("------------ up")
  print(one_comp)
  one_de <- read.csv(sprintf("out.de/%s",one_comp))
  #glist <- as.character(na.omit(one_de$X[one_de$pvalue<pcutoff & one_de$log2FoldChange>0]))  
  #enriched <- enrichr(glist, dbs)
  glist <- as.character(na.omit(one_de$X[one_de$log2FoldChange>0]))[1:1000]
  enriched <- enrichr(glist, dbs)
  
  for(cur_db in dbs){
    
    ###### normal GO plotting
    df <- enriched[[cur_db]]
    df <- df[df$P.value<0.1,]
    df <- df[order(df$P.value,decreasing = TRUE),]
    write.csv(df, sprintf("newout/up50_%s_%s.csv",one_comp, cur_db))
    df$logp <- -log10(df$P.value)
    df$cond <- one_comp
    df$cur_db <- cur_db
    newoutgo <- rbind(newoutgo, df)
    
    
    percpos <- function(x) mean(one_de$log2FoldChange[one_de$X %in% str_split(x,";")[[1]]]>0)
    df$percpos <- sapply(df$Genes,percpos)
    
    countgogenes <- function(x) length(str_split(x,";")[[1]])
    df$numgene <- sapply(df$Genes,countgogenes)
    
    df1 <- df
    df1$part_logp <- df$logp*df$percpos
    df1$isup <- "up"
    df2 <- df
    df2$part_logp <- df$logp*(1-df$percpos)
    df2$isup <- "down"
    df <- rbind(df1,df2)
    
    list_all_plot <- list()
    df <- df[df$numgene>1,]
    if(nrow(df)>0){
      df$Term <- str_split_fixed(df$Term,"\\(",2)[,1]
      df$Term <- str_split_fixed(df$Term,"_",2)[,1]
      
      df$db <- cur_db
      
      alldf <- df
      df$Term <- str_trim(df$Term)
      
      df <- df[order(df$logp,decreasing = FALSE),]
      df$Term<-factor(df$Term, levels = unique(df$Term))
      
      df <- df[df$Adjusted.P.value<1e-2,]
      
      p <- ggplot(df, aes_string(x = "Term", y = "part_logp", fill="isup")) + 
        geom_bar(stat = "identity") + coord_flip() + 
        xlab("Downregulated")+
        ylab("-Log10(p-value)")+
        theme_bw()
      p
      ggsave(sprintf("newout/go_up50_%s_%s.pdf",one_comp,cur_db), width = 10, height = 0.15*nrow(df))      
    }
    allGOdf <- rbind(allGOdf, df)
    
    
  }
  
  #table_go <- rbind(table_go, c(one_comp, "log2FoldChange>0", "", "", ""))
  
  print("------------ down")
  print(one_comp)
  one_de <- read.csv(sprintf("out.de/%s",one_comp))
  #glist <- as.character(na.omit(one_de$X[one_de$padj<pcutoff & one_de$log2FoldChange<0]))  
  #enriched <- enrichr(glist, dbs)
  glist <- as.character(na.omit(one_de$X[one_de$log2FoldChange<0]))[1:1000]
  enriched <- enrichr(glist, dbs)
  
  for(cur_db in dbs){
    
    ###### normal GO plotting
    df <- enriched[[cur_db]]
    df <- df[df$P.value<0.1,]
    df <- df[order(df$P.value,decreasing = TRUE),]
    write.csv(df, sprintf("newout/neg50_%s_%s.csv",one_comp, cur_db))
    df$logp <- -log10(df$P.value)
    df$cond <- one_comp
    df$cur_db <- cur_db
    newoutgo <- rbind(newoutgo, df)
    
    percpos <- function(x) mean(one_de$log2FoldChange[one_de$X %in% str_split(x,";")[[1]]]>0)
    df$percpos <- sapply(df$Genes,percpos)
    
    countgogenes <- function(x) length(str_split(x,";")[[1]])
    df$numgene <- sapply(df$Genes,countgogenes)
    
    df1 <- df
    df1$part_logp <- df$logp*df$percpos
    df1$isup <- "up"
    df2 <- df
    df2$part_logp <- df$logp*(1-df$percpos)
    df2$isup <- "down"
    df <- rbind(df1,df2)
    
    df <- df[df$numgene>1,]
    
    if(nrow(df)>0){
      
      df$Term <- str_split_fixed(df$Term,"\\(",2)[,1]
      df$Term <- str_split_fixed(df$Term,"_",2)[,1]
      
      df$db <- cur_db
      
      alldf <- df
      df$Term <- str_trim(df$Term)
      
      df <- df[order(df$logp,decreasing = FALSE),]
      df$Term<-factor(df$Term, levels = unique(df$Term))
      
      df <- df[df$Adjusted.P.value<1e-2,]
      
      p <- ggplot(df, aes_string(x = "Term", y = "part_logp", fill="isup")) + 
        geom_bar(stat = "identity") + coord_flip() + 
        xlab("Upregulated")+
        ylab("-Log10(p-value)")+
        theme_bw()
      p
      ggsave(sprintf("newout/go_neg50_%s_%s.pdf",one_comp,cur_db), width = 10, height = 0.15*nrow(df))
    }
    allGOdf <- rbind(allGOdf, df)
  }
}

