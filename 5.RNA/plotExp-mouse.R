#!/usr/bin/Rscript
source("~/R/library.R")
genes=commandArgs(T)
#genes="Nrf2"
gene.m=genes[1]

read.tab("~/workspace/e.Blastocyst-K27ac/5.RNA/1.Stringtie/NT-NF.reps.TPM.txt",row.names = 1,header=T) ->cdata
cdat2 = cdata[row.names(cdata)==gene.m,grep("^NF.*_rep[1-9][1-9]?$",colnames(cdata))] %>% t %>% `colnames<-`("TPM")
rownames(cdat2) %>% str.split(.,"[-_]",col.names = c("sample","stage","rep")) %>% cbind(.,cdat2) -> cdat3 
cdat3$stage %<>% factor(levels = c("oocyte","zygote","e2cell","l2cell","4cell","8cell","Morula","ICM","TE"))

temp = data.table(cdat3)
temp = temp[,.(mean=mean(TPM),sem=sd(TPM)/sqrt(length(TPM))),by=.(sample,stage)]
temp = data.frame(temp)

ggplot() +  geom_bar(data=temp,aes(x=stage,weight=mean),fill="skyblue",
                     position = position_dodge2(width = 0.5),width = 0.8,linewidth=0.75) + 
    geom_errorbar(data=temp,aes(x=stage, ymin=mean-sem, ymax=mean+sem), 
                  width=0.4,alpha=1) +
    geom_jitter(data=cdat3,aes(x=stage,y=TPM,fill=stage),
                position = position_jitterdodge(jitter.width = 0.4),linewidth=0.5) + 
    ggtitle(gene.m) + xlab("") + theme_bw() + 
    theme(aspect.ratio = 1,axis.text.x=element_text(angle = 45,hjust=1))-> g2
ggplot2::ggsave(paste0("~/workspace/e.Blastocyst-K27ac/5.RNA/1.Stringtie/5.jitter-gene-TPM/",gene.m,".mouse.pdf"),
                g2,width = 5,height = 5,device = cairo_pdf)


