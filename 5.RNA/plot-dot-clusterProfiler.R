#!/usr/bin/Rscript
source("/home/qszhu/R/library.R")
library(ggplot2)


#samples =1:7#c("Oocyte","rebuild","inherit","erase")
evalue="qvalue" #展示pvalue或evalue
sel.num = 3 #每个cluster选择其中最多3个evalue最小的term
max.value = 0.05 #每个cluster最大可接受的evalue的值
min.count = 5
max.e = 30

paste0(wdir,"/",key,".txt") %>% read.tab -> cdata
if(T){
    cdata$Cluster %<>% factor()
}else{
    cdata$Cluster %<>% factor(levels=samples)
}

cdata[,c("Description",evalue,"Cluster")] %>% reshape2::dcast(.,Description~Cluster,value.var=evalue) -> cdata2
cdata2[is.na(cdata2)] = 1

cdata[,c("Description","Count","Cluster")] %>% reshape2::dcast(.,Description~Cluster,value.var="Count") -> cdata3
cdata3[is.na(cdata3)] = 0



sel.terms = c()
for(i in 2:length(cdata2)){
    cdata2[order(cdata2[,i]),] -> temp
    temp[temp[,i]<max.value,] -> temp
    cdata3[cdata3[,i]>min.count,1] -> temp2
    temp[1:sel.num,1] -> st
    st[st %in% temp2]-> st
    sel.terms %<>% c(.,st)
}

sel.terms %<>% .[!is.na(.)] %>% .[!duplicated(.)]

cdata3[which(cdata3$Description %in% sel.terms),] %>%
    reshape2::melt(.,id.var="Description",value.name="Count",variable.name="Cluster") -> cdata5

cdata2[which(cdata2$Description %in% sel.terms),] ->  cdata4
reshape2::melt(cdata4,id.var="Description",value.name=evalue,variable.name="Cluster") -> cdata4

merge(cdata4,cdata5) -> cdata4

cdata4[,evalue] %<>% log10 %>% multiply_by(-1)
cdata4[,evalue][cdata4[,evalue]>max.e] = max.e
rang=range(cdata4[,evalue]) %>% c(5,.) %>% sort

cdata4 = cdata4[cdata4$Count!=0,]
#将count数目规整到一系列数，取到区间下沿，再用size确定他们的大小
j=-1
sizes=c(1,2,3,4,5)
for(i in c(100,140,180,220,260)){ #c(1,5,10,20,40)
    cdata4$Size[cdata4$Count>j&cdata4$Count<=i] = i
    j=i
}
cdata4$Size[cdata4$Count>j] = j
cdata4$Size %<>% as.factor()


cdata4$Description %<>% factor(levels=sel.terms[length(sel.terms):1])

if(T){
    cdata4$Cluster %<>% factor()
}else{
    cdata4$Cluster %<>% factor(levels=samples)
}

write.tab(cdata4,paste0(wdir,"/",key,".",evalue,".tab"))
ggplot(cdata4) + geom_point(aes_string(y="Description",x="Cluster",size="Size",color=evalue)) +
    scale_color_gradientn(colors = c("darkblue","grey","red"),
                          values = scales::rescale(rang,from=range(rang),to=c(0,1)),
                          name=paste0("-log10(",evalue,")")) +
    scale_size_manual(values = sizes,name="gene count") +
    scale_x_discrete(position = "top") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = -0.1),
          panel.grid.major = element_blank(),
          aspect.ratio = 2)

ggplot2::ggsave(paste0(wdir,"/",key,".",evalue,".dot.pdf"),width = 10,height = 10,useDingbats =F)





