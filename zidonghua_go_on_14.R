
##C:\1sci\Aceshi4

dir="c:\\1sci\\Aceshi4\\CESHI"

setwd(dir)

# install.packages("dplyr")

library(dplyr)
###################################################################    gene_id

library(pheatmap)
library(limma)
library(grid)
library(DESeq2)
###最原始的数据可能第一列不是gene_id，为方便下一步除重工作需要更改第一列行名为gene_id


###最原始的数据可能存在重复值，第一步清除重复值导出数据DEG    data0_0为原始数据的名称
data0_00<-read.csv("cancer_covid_Calu3.csv.",header=T)#插入待分析的文本   #$$  gene_id

data0_00 %>% 
  distinct(gene_id,.keep_all=T) ->data0_000

x0<-Reduce(intersect,list(data0_000$gene_id))

rownames(data0_000)<-data0_000$gene_id


data0_000[c(1,1000000000),]


data0_0<-data0_000[x0,-1]



                          
write.csv(data0_0,file="data0_0.csv")


Sys.sleep(10)



###################################################################
#加入各种天花乱坠的操作
###################################################################
###################################################################
###################################################################################################################### stats是什么包



setwd(dir)

# tcga<-read.table("data0_0.csv.txt",header = T,row.names = 1)     





tcga<-read.csv("data0_0.csv.",header = T,row.names = 1) #此处报错'row.names'里不能有遗漏值可能为NAMEN有NA



#raw.count<-read.table("SUM.csv",header = T,row.names = 1)
#进行条件设置
condition<-c(rep('Tumor',3),rep('Normal',3))######################对照组改数@@
coldata<- data.frame(row.names=colnames(tcga), condition)

dds <- DESeqDataSetFromMatrix(countData = tcga,
                              colData = coldata,
                              design = ~ condition)
#设置对照 
dds$condition <- relevel(dds$condition,'Normal')
#计算开始，样本量大的话，可以先去干点别的
dds <- DESeq(dds)
#get结果
res <- results(dds)
#设置cutoff，表达差异大于1或者小于-1，且padj<0.01
resSig <- subset(res, abs(log2FoldChange)>1 & padj < 0.01)
#输出结果
resSig<-data.frame(resSig)
write.csv(resSig,file="DEG.csv")

###################################################
# tcga=as.matrix(tcga)
# 
# rownames(tcga)=tcga[,1]
##此处报错可能是数值格式不统一 或者表格基因名带#  有可能是标题行有重复名称
GeneExp=tcga[,1:ncol(tcga)] #此位置不可为2
TCGA=matrix(as.numeric(as.matrix(GeneExp)),nrow=nrow(GeneExp),dimnames=list(rownames(GeneExp),colnames(GeneExp)))
TCGA=avereps(TCGA)

#过滤表达量低的基???
# TCGA=TCGA[rowMeans(TCGA)>1,]
# TCGA=round(TCGA,0)


desdsf2 <- DESeq(dds, parallel = T)
coldata<- data.frame(row.names=colnames(tcga), condition)
#此处要注意raw.count的排序需要与condition顺序一致
#构建deseq2对象


Allgene<- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
Allgene=na.omit(Allgene)

write.csv(resSig,file="DEG1.csv")
write.table(Allgene,"allgene.txt",sep="\t",quote = F)
###############################
#保存差异基因
Diffgene = Allgene[(Allgene$padj < 0.05 & (Allgene$log2FoldChange>=1 | Allgene$log2FoldChange<=(-1))),]
write.table(Diffgene, "Diffgene.txt",sep="\t",quote=F)
write.csv(resSig,file="Diffgene.csv")
#保存上调的差异基因
Upgene = Allgene[(Allgene$padj < 0.05 & (Allgene$log2FoldChange>=1)),]
write.table(Upgene, "Upgene.txt",sep="\t",quote=F)
#保存下调的差异基因
Downgene = Allgene[(Allgene$padj < 0.05 & (Allgene$log2FoldChange<=(-1))),]
write.table(Downgene, "Downgene.txt",sep="\t",quote=F)
#保存所有矫正后基因的表达量
Normalizegeneexp=as.data.frame(counts(desdsf2 , normalized=TRUE)) 
Normalizegeneexp1=rbind(id=colnames(Normalizegeneexp),Normalizegeneexp)
write.table(Normalizegeneexp1,"Normalizegeneexp.txt",sep="\t",quote=F,col.names=F)   
#保存差异基因的表达量
Diffgeneexp=rbind(id=colnames(Normalizegeneexp1),Normalizegeneexp1[rownames(Diffgene),])
write.table(Diffgeneexp,"Diffgeneexp.txt",sep="\t",quote=F,col.names=F)

# #绘制热图
# inputheatmap<-Normalizegeneexp[rownames(Diffgene),]
# inputheatmap=log2(inputheatmap+1)
# inputheatmap=inputheatmap[1:20,]
# Type=c(rep("group1",3),rep("group2",3))#############################对照组改数@@
# names(Type)=colnames(inputheatmap)
# Type=as.data.frame(Type)
# ann_colors = list(
#   group = c(group1="green", group2="red"))
# 
# pdf("heatmap.pdf",10,8)
# pheatmap(inputheatmap, annotation_col=Type,
#          color = colorRampPalette(c("green", "black", "red"))(50),cluster_cols =F,
#          fontsize = 10,fontsize_row=10,fontsize_col=5,
#          show_colnames = F,
#          annotation_legend = T,
#          annotation_names_col = T,
#          annotation_colors =ann_colors[1])

# 
# dev.off()

#绘制火山???
nosig<-Allgene[abs(Allgene$log2FoldChange)< 1 | Allgene$padj>=0.05,]
xmax<-max(Allgene$log2FoldChange)
ymax<-max(-log10(Allgene$padj))
downgene<- transform(Downgene,padj=-log10(Downgene$padj))
upgene<- transform(Upgene,padj=-log10(Upgene$padj))
nosig<- transform(nosig,padj=-log10(nosig$padj))



pdf("Volcano.pdf")
plot(nosig$log2FoldChange,nosig$padj,xlim = c(-xmax,xmax),ylim=c(0,900),col="black",
     pch=16,cex=0.9,main = "Volcano",xlab = "log2FoldChange",ylab="-log10(padj)")


points(upgene$log2FoldChange,upgene$padj,col="red",pch=16,cex=0.9)

points(downgene$log2FoldChange,downgene$padj,col="green",pch=16,cex=0.9)



abline(v=0,lwd=3,lty=2)
dev.off()


# 
###################################################################
###################################################################


###################################################################
#DEG基因名

# DEG.csv
# 
# colnames(df)[1] <- 'bb'




###################################################################
#数据导出全癌基因222222222222222222222222222222222222222222222222222222%%
###################################################################%%

#数据导出必为癌驱动基因CSV################################################
data1<-read.csv("DEG.csv",header=T)  #测试版


colnames(data1)[1] <- 'gene_id'

data3<-read.csv('A1cancer_driver_mustbe.csv',header=T)#

x2<-Reduce(intersect,list(data1$gene_id,
                          data3$gene_id))

########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

data1_mustbe<-data1[x2,-1]

# write.csv(data1_mustbe,file="+A1cancerdriver_mustbe_1.csv")
####$$$$$$$$$$

write.csv(data1_mustbe,file="XA1.csv")

dataXA1<-read.csv("XA1.csv",header=T)  #测试版

colnames(dataXA1)[1] <- 'gene_id'

write.csv(dataXA1,file="+A1cancerdriver_mustbe_1.csv",row.names = FALSE)
####$$$$$$$$$$



##
# 
# data101<-read.csv("+cancerdriver_mustbe_1.csv",header=T)  #测试版
# 
# colnames(data101)[1] <- 'gene_id'
# 
# write.csv(data101,file="-01cancerdriver_mustbe_1.csv")






data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data2<-read.csv('A2cancerdriver_all.csv',header=T)#

x1<-Reduce(intersect,list(data1$gene_id,
                          data2$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

data1_1<-data1[x1,-1]

# write.csv(data1_1,file="+A2cancerdriver_all_1_1.csv")
###

write.csv(data1_1,file="XA2.csv")

dataXA2<-read.csv("XA2.csv",header=T)  #测试版

colnames(dataXA2)[1] <- 'gene_id'

write.csv(dataXA2,file="+A2cancerdriver_all_1_1.csv",row.names = FALSE)

###

# data202<-read.csv("+cancerdriver_all_1_1.csv",header=T)  #测试版
# 
# colnames(data202)[1] <- 'gene_id'
# 
# write.csv(data202,file="-02cancerdriver_all_1_1.csv")



######
#数据导出可能为癌驱动基因csv###############################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data4<-read.csv('A3cancer_driver_maybe.csv',header=T)#

x3<-Reduce(intersect,list(data1$gene_id,
                          data4$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

data1_maybe<-data1[x3,-1]

# write.csv(data1_maybe,file="+A3cancerdriver_maybe_1.csv")
###

write.csv(data1_maybe,file="XA3.csv") ###XA3

dataXA3<-read.csv("XA3.csv",header=T)  #测试版 ###XA3

colnames(dataXA3)[1] <- 'gene_id' ###XA3

write.csv(dataXA3,file="+A3cancerdriver_maybe_1.csv",row.names = FALSE) ###XA3

###

##
# 
# data303<-read.csv("+cancerdriver_maybe_1.csv",header=T)  #测试版
# 
# colnames(data303)[1] <- 'gene_id'
# 
# write.csv(data303,file="-03cancerdriver_maybe_1.csv")

######
#数据导出细胞凋亡基因csv############################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data5<-read.csv('A4CELL_DEATH.csv',header=T)#

x4<-Reduce(intersect,list(data1$gene_id,
                          data5$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

cell_death_1<-data1[x4,-1]

# write.csv(cell_death_1.csv,file="+A4cell_death_1.csv")

###

write.csv(cell_death_1,file="XA4.csv")

dataXA4<-read.csv("XA4.csv",header=T)  #测试版

colnames(dataXA4)[1] <- 'gene_id'

write.csv(dataXA4,file="+A4cell_death_1.csv",row.names = FALSE)

###



##

# data404<-read.csv("+cell_death_1.csv",header=T)  #测试版
# 
# colnames(data404)[1] <- 'gene_id'
# 
# write.csv(data404,file="-04cell_death_1.csv")

######
#数据导出抑癌基因csv#####################################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data6<-read.csv('A5anticancer.csv',header=T)#

x5<-Reduce(intersect,list(data1$gene_id,
                          data6$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

data1_anticancer<-data1[x5,-1]

# write.csv(data1_anticancer,file="+A5canticancer_1.csv")


###

write.csv(data1_anticancer,file="XA5.csv")

dataXA5<-read.csv("XA5.csv",header=T)  #测试版

colnames(dataXA5)[1] <- 'gene_id'

write.csv(dataXA5,file="+A5canticancer_1.csv",row.names = FALSE)

###
##

# data505<-read.csv("+anticancer_1.csv",header=T)  #测试版
# 
# colnames(data505)[1] <- 'gene_id'
# 
# write.csv(data505,file="-05anticancer_1.csv")
######
#数据导出细胞分裂csv###################################################
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data7<-read.csv('A6cell_cycle.csv',header=T)#

x6<-Reduce(intersect,list(data1$gene_id,
                          data7$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

data1_cell_cycle<-data1[x6,-1]

# write.csv(data1_cell_cycle,file="+A6cell_cycle_1.csv")

###

write.csv(data1_cell_cycle,file="XA6.csv") ##2

dataXA6<-read.csv("XA6.csv",header=T)  #测试版 ##2

colnames(dataXA6)[1] <- 'gene_id' ##1

write.csv(dataXA6,file="+A6cell_cycle_1.csv",row.names = FALSE) ## 2
 
###

##

#data606<-read.csv("+A6cell_cycle_1.csv",header=T)  #测试版

# colnames(data606)[1] <- 'gene_id' 
# 
# write.csv(data606,file="-06cell_cycle_1.csv")

#############################ZHEN完结撒花#####################史诗

#############史诗############

#数据导出(命名)csv###################################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data9<-read.csv('A7generepairall.csv',header=T)#

x9<-Reduce(intersect,list(data1$gene_id,
                          data9$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A7generepairall<-data1[x9,-1]

# write.csv(A7generepairall,file="+A7generepairall.csv")

###

write.csv(A7generepairall,file="XA9.csv") ##2

dataXA9<-read.csv("XA9.csv",header=T)  #测试版 ##2

colnames(dataXA9)[1] <- 'gene_id' ##1

write.csv(dataXA9,file="+A7generepairall.csv",row.names = FALSE) ## 2

###



#数据导出(命名)csv###################################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data10<-read.csv('A7generepairright.csv',header=T)#

x10<-Reduce(intersect,list(data1$gene_id,
                          data10$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A7generepairright<-data1[x10,-1]

# write.csv(A7generepairright,file="+A7generepairright.csv")


###

write.csv(A7generepairall,file="XA10.csv") ##2

dataXA10<-read.csv("XA10.csv",header=T)  #测试版 ##2

colnames(dataXA10)[1] <- 'gene_id' ##1

write.csv(dataXA10,file="+A7generepairright.csv",row.names = FALSE) ## 2

###



#数据导出(命名)csv###################################################

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data11<-read.csv('A8causing genesall.csv',header=T)#

x11<-Reduce(intersect,list(data1$gene_id,
                           data11$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A8causinggenesall<-data1[x11,-1]

# write.csv(A8causinggenesall,file="+A8causing genesall.csv")

###

write.csv(A8causinggenesall,file="X11.csv") ##2

dataX11<-read.csv("X11.csv",header=T)  #测试版 ##2

colnames(dataX11)[1] <- 'gene_id' ##1

write.csv(dataX11,file="+A8causing genesall.csv",row.names = FALSE) ## 2

###



#数据导出(命名)csv###################################################12

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data12<-read.csv('A8causing genesB1someright.csv',header=T)##$

x12<-Reduce(intersect,list(data1$gene_id,
                           data12$gene_id))########类型的共同项目

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A8causinggenesB1someright<-data1[x12,-1]##$$

# write.csv(A8causinggenesB1someright,file="+A8causing genesB1someright.csv") ##$$


###

write.csv(A8causinggenesB1someright,file="X12.csv") ##2

dataX12<-read.csv("X12.csv",header=T)  #测试版 ##2

colnames(dataX12)[1] <- 'gene_id' ##1

write.csv(dataX12,file="+A8causing genesB1someright.csv",row.names = FALSE) ## 2

###


#数据导出(命名)csv###################################################13##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data13<-read.csv('A9DNA damage sensor genesB1all.csv',header=T)##$$DATA和文档名字要改

x13<-Reduce(intersect,list(data1$gene_id,
                           data13$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A9DNAdamagesensorgenesB1all<-data1[x13,-1]##$$文档名字要改 X后序号要改

# write.csv(A9DNAdamagesensorgenesB1all,file="+A9DNA damage sensor genesB1all.csv") ##$$两个文档名字要改


###

write.csv(A9DNAdamagesensorgenesB1all,file="X13.csv") ##2

dataX13<-read.csv("X13.csv",header=T)  #测试版 ##2

colnames(dataX13)[1] <- 'gene_id' ##1

write.csv(dataX13,file="+A9DNA damage sensor genesB1all.csv",row.names = FALSE) ## 2

###

#数据导出(命名)csv###################################################14##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data14<-read.csv('A9DNA damage sensor genesB1someright.csv',header=T)##$$DATA和文档名字要改

x14<-Reduce(intersect,list(data1$gene_id,
                           data14$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A9DNAdamagesensorgenesB1someright<-data1[x14,-1]##$$文档名字要改 X后序号要改

# write.csv(A9DNAdamagesensorgenesB1someright,file="+A9DNA damage sensor genesB1someright.csv") ##$$两个文档名字要改


###

write.csv(A9DNAdamagesensorgenesB1someright,file="X14.csv") ##2

dataX14<-read.csv("X14.csv",header=T)  #测试版 ##2

colnames(dataX14)[1] <- 'gene_id' ##1

write.csv(dataX14,file="+A9DNA damage sensor genesB1someright.csv",row.names = FALSE) ## 2

###

#数据导出(命名)csv###################################################15##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data15<-read.csv('A10DNA repair pathwaysB1all.csv',header=T)##$$DATA和文档名字要改

x15<-Reduce(intersect,list(data1$gene_id,
                           data15$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A10DNArepairpathwaysB1all<-data1[x15,-1]##$$文档名字要改 X后序号要改

# write.csv(A10DNArepairpathwaysB1all,file="+A10DNA repair pathwaysB1all.csv") ##$$两个文档名字要改


###

write.csv(A10DNArepairpathwaysB1all,file="X15.csv") ####################2

dataX15<-read.csv("X15.csv",header=T)  #测试版 #########################2

colnames(dataX15)[1] <- 'gene_id' ######################################1

write.csv(dataX15,file="+A10DNA repair pathwaysB1all.csv",row.names = FALSE) ##### 2

###




#数据导出(命名)csv###################################################16##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data16<-read.csv('A10DNA repair pathwaysB1someright.csv',header=T)##$$DATA和文档名字要改

x16<-Reduce(intersect,list(data1$gene_id,
                           data16$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A10DNArepairpathwaysB1someright<-data1[x16,-1]##$$文档名字要改 X后序号要改

# write.csv(A10DNArepairpathwaysB1someright,file="+A10DNA repair pathwaysB1someright.csv") ##$$两个文档名字要改


###

write.csv(A7generepairall,file="x16.csv") ##2

datax16<-read.csv("x16.csv",header=T)  #测试版 ##2

colnames(datax16)[1] <- 'gene_id' ##1

write.csv(datax16,file="+A10DNA repair pathwaysB1someright.csv",row.names = FALSE) ## 2

###


#数据导出(命名)csv###################################################17##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data17<-read.csv('A11DNA damage response pathway.B1someright.csv',header=T)##$$DATA和文档名字要改

x17<-Reduce(intersect,list(data1$gene_id,
                           data17$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A11<-data1[x17,-1]##$$文档名字要改 X后序号要改

# write.csv(A11,file="+A11DNA damage response pathway.B1someright.csv") ##$$两个文档名字要改


###

write.csv(A11,file="x17.csv") ##2

datax17<-read.csv("x17.csv",header=T)  #测试版 ##2

colnames(datax17)[1] <- 'gene_id' ##1

write.csv(datax17,file="+A11DNA damage response pathway.B1someright.csv",row.names = FALSE) ## 2

###


#数据导出(命名)csv###################################################18##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data18<-read.csv('A11DNA damage response pathwayB1all.csv',header=T)##$$DATA和文档名字要改

x18<-Reduce(intersect,list(data1$gene_id,
                           data18$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A18<-data1[x18,-1]##$$文档名字要改 X后序号要改

# write.csv(A18,file="+A11DNA damage response pathwayB1all.csv") ##$$两个文档名字要改


###

write.csv(A18,file="x18.csv") ##2

datax18<-read.csv("x18.csv",header=T)  #测试版 ##2

colnames(datax18)[1] <- 'gene_id' ##1

write.csv(datax18,file="+A11DNA damage response pathwayB1all.csv",row.names = FALSE) ## 2

###



#数据导出(命名)csv###################################################19##$

data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data19<-read.csv('A12ATR-ATM pathway.B1all.csv',header=T)##$$DATA和文档名字要改

x19<-Reduce(intersect,list(data1$gene_id,
                           data19$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A19<-data1[x19,-1]##$$文档名字要改 X后序号要改

# write.csv(A19,file="+A12ATR-ATM pathway.B1all.csv") ##$$两个文档名字要改


###

write.csv(A19,file="x19.csv") ##2

datax19<-read.csv("x19.csv",header=T)  #测试版 ##2

colnames(datax19)[1] <- 'gene_id' ##1

write.csv(datax19,file="+A12ATR-ATM pathway.B1all.csv",row.names = FALSE) ## 2

###


#数据导出(命名)csv###################################################20##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data20<-read.csv('A12ATR-ATM pathway.B1someright.csv',header=T)##$$DATA和文档名字要改

x20<-Reduce(intersect,list(data1$gene_id,
                           data20$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A20<-data1[x20,-1]##$$文档名字要改 X后序号要改

# write.csv(A20,file="+A12ATR-ATM pathway.B1someright.csv") ##$$两个文档名字要改


###

write.csv(A20,file="x20.csv") ##2

datax20<-read.csv("x20.csv",header=T)  #测试版 ##2

colnames(datax20)[1] <- 'gene_id' ##1

write.csv(datax20,file="+A12ATR-ATM pathway.B1someright.csv",row.names = FALSE) ## 2

###


#数据导出(命名)csv###################################################21##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data21<-read.csv('A13DNA damage checkpoint.B1all.csv',header=T)##$$DATA和文档名字要改

x21<-Reduce(intersect,list(data1$gene_id,
                           data21$gene_id))########类型的共同项目##$$两个X后需要要改

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A21<-data1[x21,-1]##$$文档名字要改 X后序号要改

# write.csv(A21,file="+A13DNA damage checkpoint.B1all.csv") ##$$两个文档名字要改


###

write.csv(A21,file="x21.csv") ##2

datax21<-read.csv("x21.csv",header=T)  #测试版 ##2

colnames(datax21)[1] <- 'gene_id' ##1

write.csv(datax21,file="+A13DNA damage checkpoint.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################22##$  ???
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data22<-read.csv('A13DNA damage checkpoint.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM

x22<-Reduce(intersect,list(data1$gene_id,
                           data22$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A22<-data1[x22,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A22,file="+A13DNA damage checkpoint.B1someright.csv") ##$$两个文档名字要改SSSSSSSMMMMMMM


###

write.csv(A22,file="x22.csv") ##2

datax22<-read.csv("x22.csv",header=T)  #测试版 ##2

colnames(datax22)[1] <- 'gene_id' ##1

write.csv(datax22,file="+A13DNA damage checkpoint.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################23##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data23<-read.csv('A14DNA damage signaling.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM

x23<-Reduce(intersect,list(data1$gene_id,
                           data23$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A23<-data1[x23,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A23,file="+A14DNA damage signaling.B1all.csv") ##$$两个文档名字要改SSSSSSSMMMMMMM


###

write.csv(A23,file="x23.csv") ##2

datax23<-read.csv("x23.csv",header=T)  #测试版 ##2

colnames(datax23)[1] <- 'gene_id' ##1

write.csv(datax23,file="+A14DNA damage signaling.B1all.csv",row.names = FALSE) ## 2

###



##################数据导出(命名)csv###################################################24##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data24<-read.csv('A14DNA damage signaling.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM

x24<-Reduce(intersect,list(data1$gene_id,
                           data24$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A24<-data1[x24,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

write.csv(A24,file="+A14DNA damage signaling.B1someright.csv") ##$$两个文档名字要改SSSSSSSMMMMMMM


###

write.csv(A24,file="x24.csv") ##2

datax24<-read.csv("x24.csv",header=T)  #测试版 ##2

colnames(datax24)[1] <- 'gene_id' ##1

write.csv(datax24,file="+A14DNA damage signaling.B1someright.csv",row.names = FALSE) ## 2

###



##################数据导出(命名)csv###################################################25##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data25<-read.csv('A15DNA damage response (DDR).B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM

x25<-Reduce(intersect,list(data1$gene_id,
                           data25$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A25<-data1[x25,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A25,file="+A15DNA damage response (DDR).B1all.csv") ##$$两个文档名字要改SSSSSSSMMMMMMM


###

write.csv(A25,file="x25.csv") ##2

datax25<-read.csv("x25.csv",header=T)  #测试版 ##2

colnames(datax25)[1] <- 'gene_id' ##1

write.csv(datax25,file="+A15DNA damage response (DDR).B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################26##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data26<-read.csv('A15DNA damage response (DDR).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x26<-Reduce(intersect,list(data1$gene_id,
                           data26$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A26<-data1[x26,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A26,file="+A15DNA damage response (DDR).B1someright.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A26,file="x26.csv") ##2

datax26<-read.csv("x26.csv",header=T)  #测试版 ##2

colnames(datax26)[1] <- 'gene_id' ##1

write.csv(datax26,file="+A15DNA damage response (DDR).B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################27##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data27<-read.csv('A17telomere shortening.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x27<-Reduce(intersect,list(data1$gene_id,
                           data27$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A27<-data1[x27,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A27,file="+A17telomere shortening.B1all.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A27,file="x27.csv") ##2

datax27<-read.csv("x27.csv",header=T)  #测试版 ##2

colnames(datax27)[1] <- 'gene_id' ##1

write.csv(datax27,file="+A17telomere shortening.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################28##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data28<-read.csv('A17telomere shortening.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x28<-Reduce(intersect,list(data1$gene_id,
                           data28$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A28<-data1[x28,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

write.csv(A28,file="+A17telomere shortening.B1someright.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A28,file="x28.csv") ##2

datax28<-read.csv("x28.csv",header=T)  #测试版 ##2

colnames(datax28)[1] <- 'gene_id' ##1

write.csv(datax28,file="+A17telomere shortening.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################29##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data29<-read.csv('A18telomere length regulation genes.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x29<-Reduce(intersect,list(data1$gene_id,
                           data29$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A29<-data1[x29,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

write.csv(A29,file="+A18telomere length regulation genes.B1all.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1

###

write.csv(A29,file="x29.csv") ##2

datax29<-read.csv("x29.csv",header=T)  #测试版 ##2

colnames(datax29)[1] <- 'gene_id' ##1

write.csv(datax29,file="+A18telomere length regulation genes.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################30##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data30<-read.csv('A18telomere length regulation genes.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x30<-Reduce(intersect,list(data1$gene_id,
                           data30$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A30<-data1[x30,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A30,file="+A18telomere length regulation genes.B1someright.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A30,file="x30.csv") ##2

datax30<-read.csv("x30.csv",header=T)  #测试版 ##2

colnames(datax30)[1] <- 'gene_id' ##1

write.csv(datax30,file="+A18telomere length regulation genes.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################31##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data31<-read.csv('A19telomere_associated disease genes.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x31<-Reduce(intersect,list(data1$gene_id,
                           data31$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A31<-data1[x31,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A31,file="+A19telomere_associated disease genes.B1all.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A31,file="x31.csv") ##2

datax31<-read.csv("x31.csv",header=T)  #测试版 ##2

colnames(datax31)[1] <- 'gene_id' ##1

write.csv(datax31,file="+A19telomere_associated disease genes.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################32##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data32<-read.csv('A19telomere_associated disease genes.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x32<-Reduce(intersect,list(data1$gene_id,
                           data32$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A32<-data1[x32,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A32,file="+A19telomere_associated disease genes.B1someright.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A32,file="x32.csv") ##2

datax32<-read.csv("x32.csv",header=T)  #测试版 ##2

colnames(datax32)[1] <- 'gene_id' ##1

write.csv(datax32,file="+A19telomere_associated disease genes.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################32##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data32<-read.csv('A20telomere protective protein genes.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x32<-Reduce(intersect,list(data1$gene_id,
                           data32$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A32<-data1[x32,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

write.csv(A32,file="+A20telomere protective protein genes.B1all.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A32,file="x32.csv") ##2

datax32<-read.csv("x32.csv",header=T)  #测试版 ##2

colnames(datax32)[1] <- 'gene_id' ##1

write.csv(datax32,file="+A20telomere protective protein genes.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################33##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data33<-read.csv('A20telomere protective protein genes.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x33<-Reduce(intersect,list(data1$gene_id,
                           data33$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A33<-data1[x33,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

# write.csv(A33,file="+A20telomere protective protein genes.B1all.csv") ##$$两个文档名字要改SSSSSSS 1 MMMMMMM 1


###

write.csv(A33,file="x33.csv") ##2

datax33<-read.csv("x33.csv",header=T)  #测试版 ##2

colnames(datax33)[1] <- 'gene_id' ##1

write.csv(datax33,file="+A20telomere protective protein genes.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################34##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data34<-read.csv('A21Collagen.all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x34<-Reduce(intersect,list(data1$gene_id,
                           data34$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A34<-data1[x34,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A34,file="x34.csv") ##2

datax34<-read.csv("x34.csv",header=T)  #测试版 ##2

colnames(datax34)[1] <- 'gene_id' ##1

write.csv(datax34,file="+A21Collagen.all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################35##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data35<-read.csv('A21Collagen.someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x35<-Reduce(intersect,list(data1$gene_id,
                           data35$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A35<-data1[x35,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A35,file="x35.csv") ##2

datax35<-read.csv("x35.csv",header=T)  #测试版 ##2

colnames(datax35)[1] <- 'gene_id' ##1

write.csv(datax35,file="+A21Collagen.someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################36##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data36<-read.csv('A22Fibroblast growth factor.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x36<-Reduce(intersect,list(data1$gene_id,
                           data36$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A36<-data1[x36,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A36,file="x36.csv") ##2

datax36<-read.csv("x36.csv",header=T)  #测试版 ##2

colnames(datax36)[1] <- 'gene_id' ##1

write.csv(datax36,file="+A22Fibroblast growth factor.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################37##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data37<-read.csv('A22Fibroblast growth factor.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x37<-Reduce(intersect,list(data1$gene_id,
                           data37$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A37<-data1[x37,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A37,file="x37.csv") ##2

datax37<-read.csv("x37.csv",header=T)  #测试版 ##2

colnames(datax37)[1] <- 'gene_id' ##1

write.csv(datax37,file="+A22Fibroblast growth factor.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################38##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data38<-read.csv('A23Hemoglobin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x38<-Reduce(intersect,list(data1$gene_id,
                           data38$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A38<-data1[x38,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A38,file="x38.csv") ##2

datax38<-read.csv("x38.csv",header=T)  #测试版 ##2

colnames(datax38)[1] <- 'gene_id' ##1

write.csv(datax38,file="+A23Hemoglobin.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################39##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data39<-read.csv('A23Hemoglobin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x39<-Reduce(intersect,list(data1$gene_id,
                           data39$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A39<-data1[x39,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A39,file="x39.csv") ##2

datax39<-read.csv("x39.csv",header=T)  #测试版 ##2

colnames(datax39)[1] <- 'gene_id' ##1

write.csv(datax39,file="+A23Hemoglobin.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################40##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data40<-read.csv('A24Anti-enzyme.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x40<-Reduce(intersect,list(data1$gene_id,
                           data40$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A40<-data1[x40,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A40,file="x40.csv") ##2

datax40<-read.csv("x40.csv",header=T)  #测试版 ##2

colnames(datax40)[1] <- 'gene_id' ##1

write.csv(datax40,file="+A24Anti-enzyme.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################41##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data41<-read.csv('A24Anti-enzyme.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x41<-Reduce(intersect,list(data1$gene_id,
                           data41$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A41<-data1[x41,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A41,file="x41.csv") ##2

datax41<-read.csv("x41.csv",header=T)  #测试版 ##2

colnames(datax41)[1] <- 'gene_id' ##1

write.csv(datax41,file="+A24Anti-enzyme.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################42##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data42<-read.csv('A25Kinase.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x42<-Reduce(intersect,list(data1$gene_id,
                           data42$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A42<-data1[x42,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A42,file="x42.csv") ##2

datax42<-read.csv("x42.csv",header=T)  #测试版 ##2

colnames(datax42)[1] <- 'gene_id' ##1

write.csv(datax42,file="+A25Kinase.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################43##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data43<-read.csv('A25Kinase.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x43<-Reduce(intersect,list(data1$gene_id,
                           data43$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A43<-data1[x43,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A43,file="x43.csv") ##2

datax43<-read.csv("x43.csv",header=T)  #测试版 ##2

colnames(datax43)[1] <- 'gene_id' ##1

write.csv(datax43,file="+A25Kinase.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################44##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data44<-read.csv('A26Antibodies.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x44<-Reduce(intersect,list(data1$gene_id,
                           data44$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A44<-data1[x44,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A44,file="x44.csv") ##2

datax44<-read.csv("x44.csv",header=T)  #测试版 ##2

colnames(datax44)[1] <- 'gene_id' ##1

write.csv(datax44,file="+A26Antibodies.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################45##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data45<-read.csv('A26Antibodies.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x45<-Reduce(intersect,list(data1$gene_id,
                           data45$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A45<-data1[x45,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A45,file="x45.csv") ##2

datax45<-read.csv("x45.csv",header=T)  #测试版 ##2

colnames(datax45)[1] <- 'gene_id' ##1

write.csv(datax45,file="+A26Antibodies.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################46##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data46<-read.csv('A27Hormones.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x46<-Reduce(intersect,list(data1$gene_id,
                           data46$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A46<-data1[x46,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A46,file="x46.csv") ##2

datax46<-read.csv("x46.csv",header=T)  #测试版 ##2

colnames(datax46)[1] <- 'gene_id' ##1

write.csv(datax46,file="+A27Hormones.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################47##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data47<-read.csv('A27Hormones.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x47<-Reduce(intersect,list(data1$gene_id,
                           data47$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A47<-data1[x47,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A47,file="x47.csv") ##2

datax47<-read.csv("x47.csv",header=T)  #测试版 ##2

colnames(datax47)[1] <- 'gene_id' ##1

write.csv(datax47,file="+A27Hormones.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################48##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data48<-read.csv('A28Myostatin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x48<-Reduce(intersect,list(data1$gene_id,
                           data48$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A48<-data1[x48,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A48,file="x48.csv") ##2

datax48<-read.csv("x48.csv",header=T)  #测试版 ##2

colnames(datax48)[1] <- 'gene_id' ##1

write.csv(datax48,file="+A28Myostatin.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################49##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data49<-read.csv('A28Myostatin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x49<-Reduce(intersect,list(data1$gene_id,
                           data49$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A49<-data1[x49,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A49,file="x49.csv") ##2

datax49<-read.csv("x49.csv",header=T)  #测试版 ##2

colnames(datax49)[1] <- 'gene_id' ##1

write.csv(datax49,file="+A28Myostatin.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################50##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data50<-read.csv('A29Lysozyme.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x50<-Reduce(intersect,list(data1$gene_id,
                           data50$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A50<-data1[x50,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A50,file="x50.csv") ##2

datax50<-read.csv("x50.csv",header=T)  #测试版 ##2

colnames(datax50)[1] <- 'gene_id' ##1

write.csv(datax50,file="+A29Lysozyme.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################51##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data51<-read.csv('A29Lysozyme.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x51<-Reduce(intersect,list(data1$gene_id,
                           data51$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A51<-data1[x51,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A51,file="x51.csv") ##2

datax51<-read.csv("x51.csv",header=T)  #测试版 ##2

colnames(datax51)[1] <- 'gene_id' ##1

write.csv(datax51,file="+A29Lysozyme.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################52##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data52<-read.csv('A30Immunoglobulins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x52<-Reduce(intersect,list(data1$gene_id,
                           data52$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A52<-data1[x52,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A52,file="x52.csv") ##2

datax52<-read.csv("x52.csv",header=T)  #测试版 ##2

colnames(datax52)[1] <- 'gene_id' ##1

write.csv(datax52,file="+A30Immunoglobulins.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################53##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data53<-read.csv('A30Immunoglobulins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x53<-Reduce(intersect,list(data1$gene_id,
                           data53$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A53<-data1[x53,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A53,file="x53.csv") ##2

datax53<-read.csv("x53.csv",header=T)  #测试版 ##2

colnames(datax53)[1] <- 'gene_id' ##1

write.csv(datax53,file="+A30Immunoglobulins.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################54##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data54<-read.csv('A31Extracellular matrix (ECM) proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x54<-Reduce(intersect,list(data1$gene_id,
                           data54$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A54<-data1[x54,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A54,file="x54.csv") ##2

datax54<-read.csv("x54.csv",header=T)  #测试版 ##2

colnames(datax54)[1] <- 'gene_id' ##1

write.csv(datax54,file="+A31Extracellular matrix (ECM) proteins.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################55##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data55<-read.csv('A32Insulin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x55<-Reduce(intersect,list(data1$gene_id,
                           data55$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A55<-data1[x55,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A55,file="x55.csv") ##2

datax55<-read.csv("x55.csv",header=T)  #测试版 ##2

colnames(datax55)[1] <- 'gene_id' ##1

write.csv(datax55,file="+A32Insulin.B1all.csv",row.names = FALSE) ## 2

###
####
##################数据导出(命名)csv###################################################56##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data56<-read.csv('A32Insulin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x56<-Reduce(intersect,list(data1$gene_id,
                           data56$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A56<-data1[x56,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A56,file="x56.csv") ##2

datax56<-read.csv("x56.csv",header=T)  #测试版 ##2

colnames(datax56)[1] <- 'gene_id' ##1

write.csv(datax56,file="+A32Insulin.B1someright.csv",row.names = FALSE) ## 2

###
####
##################数据导出(命名)csv###################################################57##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data57<-read.csv('A33Cholesterol transport protein (HDL_LDL).B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x57<-Reduce(intersect,list(data1$gene_id,
                           data57$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A57<-data1[x57,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A57,file="x57.csv") ##2

datax57<-read.csv("x57.csv",header=T)  #测试版 ##2

colnames(datax57)[1] <- 'gene_id' ##1

write.csv(datax57,file="+A33Cholesterol transport protein (HDL_LDL).B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################58##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data58<-read.csv('A33Cholesterol transport protein (HDL_LDL).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x58<-Reduce(intersect,list(data1$gene_id,
                           data58$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A58<-data1[x58,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A58,file="x58.csv") ##2

datax58<-read.csv("x58.csv",header=T)  #测试版 ##2

colnames(datax58)[1] <- 'gene_id' ##1

write.csv(datax58,file="+A33Cholesterol transport protein (HDL_LDL).B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################59##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data59<-read.csv('A34Protease.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x59<-Reduce(intersect,list(data1$gene_id,
                           data59$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A59<-data1[x59,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A59,file="x59.csv") ##2

datax59<-read.csv("x59.csv",header=T)  #测试版 ##2

colnames(datax59)[1] <- 'gene_id' ##1

write.csv(datax59,file="+A34Protease.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################60##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data60<-read.csv('A34Protease.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x60<-Reduce(intersect,list(data1$gene_id,
                           data60$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A60<-data1[x60,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A60,file="x60.csv") ##2

datax60<-read.csv("x60.csv",header=T)  #测试版 ##2

colnames(datax60)[1] <- 'gene_id' ##1

write.csv(datax60,file="+A34Protease.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################61##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data61<-read.csv('A35Transferrin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x61<-Reduce(intersect,list(data1$gene_id,
                           data61$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A61<-data1[x61,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A61,file="x61.csv") ##2

datax61<-read.csv("x61.csv",header=T)  #测试版 ##2

colnames(datax61)[1] <- 'gene_id' ##1

write.csv(datax61,file="+A35Transferrin.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################62##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data62<-read.csv('A35Transferrin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x62<-Reduce(intersect,list(data1$gene_id,
                           data62$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A62<-data1[x62,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A62,file="x62.csv") ##2

datax62<-read.csv("x62.csv",header=T)  #测试版 ##2

colnames(datax62)[1] <- 'gene_id' ##1

write.csv(datax62,file="+A35Transferrin.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################63##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data63<-read.csv('A36Cell membrane proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x63<-Reduce(intersect,list(data1$gene_id,
                           data63$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A63<-data1[x63,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A63,file="x63.csv") ##2

datax63<-read.csv("x63.csv",header=T)  #测试版 ##2

colnames(datax63)[1] <- 'gene_id' ##1

write.csv(datax63,file="+A36Cell membrane proteins.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################64##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data64<-read.csv('A36Cell membrane proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x64<-Reduce(intersect,list(data1$gene_id,
                           data64$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A64<-data1[x64,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A64,file="x64.csv") ##2

datax64<-read.csv("x64.csv",header=T)  #测试版 ##2

colnames(datax64)[1] <- 'gene_id' ##1

write.csv(datax64,file="+A36Cell membrane proteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################65##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data65<-read.csv('A37Receptor proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x65<-Reduce(intersect,list(data1$gene_id,
                           data65$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A65<-data1[x65,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A65,file="x65.csv") ##2

datax65<-read.csv("x65.csv",header=T)  #测试版 ##2

colnames(datax65)[1] <- 'gene_id' ##1

write.csv(datax65,file="+A37Receptor proteins.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################66##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data66<-read.csv('A37Receptor proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x66<-Reduce(intersect,list(data1$gene_id,
                           data66$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A66<-data1[x66,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A66,file="x66.csv") ##2

datax66<-read.csv("x66.csv",header=T)  #测试版 ##2

colnames(datax66)[1] <- 'gene_id' ##1

write.csv(datax66,file="+A37Receptor proteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################67##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data67<-read.csv('A38Transcription factors.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x67<-Reduce(intersect,list(data1$gene_id,
                           data67$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A67<-data1[x67,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A67,file="x67.csv") ##2

datax67<-read.csv("x67.csv",header=T)  #测试版 ##2

colnames(datax67)[1] <- 'gene_id' ##1

write.csv(datax67,file="+A38Transcription factors.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################68##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data68<-read.csv('A38Transcription factors.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x68<-Reduce(intersect,list(data1$gene_id,
                           data68$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A68<-data1[x68,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A68,file="x68.csv") ##2

datax68<-read.csv("x68.csv",header=T)  #测试版 ##2

colnames(datax68)[1] <- 'gene_id' ##1

write.csv(datax68,file="+A38Transcription factors.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################69##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data69<-read.csv('A39Protein phosphatase.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x69<-Reduce(intersect,list(data1$gene_id,
                           data69$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A69<-data1[x69,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A69,file="x69.csv") ##2

datax69<-read.csv("x69.csv",header=T)  #测试版 ##2

colnames(datax69)[1] <- 'gene_id' ##1

write.csv(datax69,file="+A39Protein phosphatase.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################70##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data70<-read.csv('A39Protein phosphatase.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x70<-Reduce(intersect,list(data1$gene_id,
                           data70$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A70<-data1[x70,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A70,file="x70.csv") ##2

datax70<-read.csv("x70.csv",header=T)  #测试版 ##2

colnames(datax70)[1] <- 'gene_id' ##1

write.csv(datax70,file="+A39Protein phosphatase.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################71##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data71<-read.csv('A40Mitochondrial proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x71<-Reduce(intersect,list(data1$gene_id,
                           data71$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A71<-data1[x71,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A71,file="x71.csv") ##2

datax71<-read.csv("x71.csv",header=T)  #测试版 ##2

colnames(datax71)[1] <- 'gene_id' ##1

write.csv(datax71,file="+A40Mitochondrial proteins.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################72##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data72<-read.csv('A40Mitochondrial proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x72<-Reduce(intersect,list(data1$gene_id,
                           data72$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A72<-data1[x72,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A72,file="x72.csv") ##2

datax72<-read.csv("x72.csv",header=T)  #测试版 ##2

colnames(datax72)[1] <- 'gene_id' ##1

write.csv(datax72,file="+A40Mitochondrial proteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################73##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data73<-read.csv('A41Oxygen oxidases.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x73<-Reduce(intersect,list(data1$gene_id,
                           data73$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A73<-data1[x73,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A73,file="x73.csv") ##2

datax73<-read.csv("x73.csv",header=T)  #测试版 ##2

colnames(datax73)[1] <- 'gene_id' ##1

write.csv(datax73,file="+A41Oxygen oxidases.B1someright.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################74##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data74<-read.csv('A42Transporter proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x74<-Reduce(intersect,list(data1$gene_id,
                           data74$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A74<-data1[x74,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A74,file="x74.csv") ##2

datax74<-read.csv("x74.csv",header=T)  #测试版 ##2

colnames(datax74)[1] <- 'gene_id' ##1

write.csv(datax74,file="+A42Transporter proteins.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################74##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data74a4<-read.csv('A42Transporter proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x74a4<-Reduce(intersect,list(data1$gene_id,
                           data74a4$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A74a4<-data1[x74a4,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A74a4,file="x74a4.csv") ##2

datax74a4<-read.csv("x74a4.csv",header=T)  #测试版 ##2

colnames(datax74a4)[1] <- 'gene_id' ##1

write.csv(datax74a4,file="+A42Transporter proteins.B1someright.csv",row.names = FALSE) ## 2
##################数据导出(命名)csv###################################################75##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data75<-read.csv('A43Cell wall proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x75<-Reduce(intersect,list(data1$gene_id,
                           data75$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A75<-data1[x75,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A75,file="x75.csv") ##2

datax75<-read.csv("x75.csv",header=T)  #测试版 ##2

colnames(datax75)[1] <- 'gene_id' ##1

write.csv(datax75,file="+A43Cell wall proteins.B1all.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################76##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data76<-read.csv('A43Cell wall proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x76<-Reduce(intersect,list(data1$gene_id,
                           data76$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A76<-data1[x76,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A76,file="x76.csv") ##2

datax76<-read.csv("x76.csv",header=T)  #测试版 ##2

colnames(datax76)[1] <- 'gene_id' ##1

write.csv(datax76,file="+A43Cell wall proteins.B1someright.csv",row.names = FALSE) ## 2

###


##################数据导出(命名)csv###################################################77##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data77<-read.csv('A44Protein phosphatase.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x77<-Reduce(intersect,list(data1$gene_id,
                           data77$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A77<-data1[x77,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A77,file="x77.csv") ##2

datax77<-read.csv("x77.csv",header=T)  #测试版 ##2

colnames(datax77)[1] <- 'gene_id' ##1

write.csv(datax77,file="+A44Protein phosphatase.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################78##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data78<-read.csv('A44Protein phosphatase.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x78<-Reduce(intersect,list(data1$gene_id,
                           data78$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A78<-data1[x78,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A78,file="x78.csv") ##2

datax78<-read.csv("x78.csv",header=T)  #测试版 ##2

colnames(datax78)[1] <- 'gene_id' ##1

write.csv(datax78,file="+A44Protein phosphatase.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################79##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data79<-read.csv('A45Nucleoproteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x79<-Reduce(intersect,list(data1$gene_id,
                           data79$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A79<-data1[x79,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A79,file="x79.csv") ##2

datax79<-read.csv("x79.csv",header=T)  #测试版 ##2

colnames(datax79)[1] <- 'gene_id' ##1

write.csv(datax79,file="+A45Nucleoproteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################80##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data80<-read.csv('A46Nucleases.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x80<-Reduce(intersect,list(data1$gene_id,
                           data80$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A80<-data1[x80,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A80,file="x80.csv") ##2

datax80<-read.csv("x80.csv",header=T)  #测试版 ##2

colnames(datax80)[1] <- 'gene_id' ##1

write.csv(datax80,file="+A46Nucleases.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################81##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data81<-read.csv('A47Hydrolases.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x81<-Reduce(intersect,list(data1$gene_id,
                           data81$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A81<-data1[x81,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A81,file="x81.csv") ##2

datax81<-read.csv("x81.csv",header=T)  #测试版 ##2

colnames(datax81)[1] <- 'gene_id' ##1

write.csv(datax81,file="+A47Hydrolases.B1all.csv",row.names = FALSE) ## 2

###
##################数据导出(命名)csv###################################################82##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data82<-read.csv('A47Hydrolases.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x82<-Reduce(intersect,list(data1$gene_id,
                           data82$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A82<-data1[x82,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A82,file="x82.csv") ##2

datax82<-read.csv("x82.csv",header=T)  #测试版 ##2

colnames(datax82)[1] <- 'gene_id' ##1

write.csv(datax82,file="+A47Hydrolases.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################83##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data83<-read.csv('A48lysosomal proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x83<-Reduce(intersect,list(data1$gene_id,
                           data83$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A83<-data1[x83,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A83,file="x83.csv") ##2

datax83<-read.csv("x83.csv",header=T)  #测试版 ##2

colnames(datax83)[1] <- 'gene_id' ##1

write.csv(datax83,file="+A48lysosomal proteins.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################84##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data84<-read.csv('A48lysosomal proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x84<-Reduce(intersect,list(data1$gene_id,
                           data84$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A84<-data1[x84,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A84,file="x84.csv") ##2

datax84<-read.csv("x84.csv",header=T)  #测试版 ##2

colnames(datax84)[1] <- 'gene_id' ##1

write.csv(datax84,file="+A48lysosomal proteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################85##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data85<-read.csv('A49extra-mitochondrial matrix proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x85<-Reduce(intersect,list(data1$gene_id,
                           data85$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A85<-data1[x85,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A85,file="x85.csv") ##2

datax85<-read.csv("x85.csv",header=T)  #测试版 ##2

colnames(datax85)[1] <- 'gene_id' ##1

write.csv(datax85,file="+A49extra-mitochondrial matrix proteins.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################86##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data86<-read.csv('A49extra-mitochondrial matrix proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x86<-Reduce(intersect,list(data1$gene_id,
                           data86$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A86<-data1[x86,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A86,file="x86.csv") ##2

datax86<-read.csv("x86.csv",header=T)  #测试版 ##2

colnames(datax86)[1] <- 'gene_id' ##1

write.csv(datax86,file="+A49extra-mitochondrial matrix proteins.B1someright.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################87##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data87<-read.csv('A50extracellular matrix proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x87<-Reduce(intersect,list(data1$gene_id,
                           data87$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A87<-data1[x87,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A87,file="x87.csv") ##2

datax87<-read.csv("x87.csv",header=T)  #测试版 ##2

colnames(datax87)[1] <- 'gene_id' ##1

write.csv(datax87,file="+A50extracellular matrix proteins.B1all.csv",row.names = FALSE) ## 2

###

##################数据导出(命名)csv###################################################88##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data88<-read.csv('A50extracellular matrix proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x88<-Reduce(intersect,list(data1$gene_id,
                           data88$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A88<-data1[x88,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A88,file="x88.csv") ##2

datax88<-read.csv("x88.csv",header=T)  #测试版 ##2

colnames(datax88)[1] <- 'gene_id' ##1

write.csv(datax88,file="+A50extracellular matrix proteins.B1someright.csv",row.names = FALSE) ## 2

###

######数据导出(命名)csv###################################################89##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data89<-read.csv('A51Antithrombin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x89<-Reduce(intersect,list(data1$gene_id,
                           data89$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A89<-data1[x89,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A89,file="x89.csv") ##2

datax89<-read.csv("x89.csv",header=T)  #测试版 ##2

colnames(datax89)[1] <- 'gene_id' ##1

write.csv(datax89,file="+A51Antithrombin.B1someright.csv",row.names = FALSE) ## 2

###
######数据导出(命名)csv###################################################89##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data89a9<-read.csv('A51Antithrombin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x89a9<-Reduce(intersect,list(data1$gene_id,
                           data89a9$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A89a9<-data1[x89a9,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A89,file="x89a9.csv") ##2

datax89a9<-read.csv("x89a9.csv",header=T)  #测试版 ##2

colnames(datax89a9)[1] <- 'gene_id' ##1

write.csv(datax89a9,file="+A51Antithrombin.B1all.csv",row.names = FALSE) ## 2
###数据导出(命名)csv###################################################90##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data90<-read.csv('A52Anti-platelet proteins.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x90<-Reduce(intersect,list(data1$gene_id,
                           data90$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A90<-data1[x90,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A90,file="x90.csv") ##2

datax90<-read.csv("x90.csv",header=T)  #测试版 ##2

colnames(datax90)[1] <- 'gene_id' ##1

write.csv(datax90,file="+A52Anti-platelet proteins.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################91##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data91<-read.csv('A52Anti-platelet proteins.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x91<-Reduce(intersect,list(data1$gene_id,
                           data91$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A91<-data1[x91,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A91,file="x91.csv") ##2

datax91<-read.csv("x91.csv",header=T)  #测试版 ##2

colnames(datax91)[1] <- 'gene_id' ##1

write.csv(datax91,file="+A52Anti-platelet proteins.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################92##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data92<-read.csv('A53Antithrombin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x92<-Reduce(intersect,list(data1$gene_id,
                           data92$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A92<-data1[x92,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A92,file="x92.csv") ##2

datax92<-read.csv("x92.csv",header=T)  #测试版 ##2

colnames(datax92)[1] <- 'gene_id' ##1

write.csv(datax92,file="+A53Antithrombin.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################93##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data93<-read.csv('A53Antithrombin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x93<-Reduce(intersect,list(data1$gene_id,
                           data93$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A93<-data1[x93,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A93,file="x93.csv") ##2

datax93<-read.csv("x93.csv",header=T)  #测试版 ##2

colnames(datax93)[1] <- 'gene_id' ##1

write.csv(datax93,file="+A53Antithrombin.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################94##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data94<-read.csv('A54Anti-thrombotic antibodies.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x94<-Reduce(intersect,list(data1$gene_id,
                           data94$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A94<-data1[x94,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A94,file="x94.csv") ##2

datax94<-read.csv("x94.csv",header=T)  #测试版 ##2

colnames(datax94)[1] <- 'gene_id' ##1

write.csv(datax94,file="+A54Anti-thrombotic antibodies.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################95##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data95<-read.csv('A55Myosin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x95<-Reduce(intersect,list(data1$gene_id,
                           data95$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A95<-data1[x95,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A95,file="x95.csv") ##2

datax95<-read.csv("x95.csv",header=T)  #测试版 ##2

colnames(datax95)[1] <- 'gene_id' ##1

write.csv(datax95,file="+A55Myosin.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################96##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data96<-read.csv('A55Myosin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x96<-Reduce(intersect,list(data1$gene_id,
                           data96$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A96<-data1[x96,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A96,file="x96.csv") ##2

datax96<-read.csv("x96.csv",header=T)  #测试版 ##2

colnames(datax96)[1] <- 'gene_id' ##1

write.csv(datax96,file="+A55Myosin.B1someright.csv",row.names = FALSE) ## 2


###数据导出(命名)csv###################################################97##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data97<-read.csv('A56Actin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x97<-Reduce(intersect,list(data1$gene_id,
                           data97$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A97<-data1[x97,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A97,file="x97.csv") ##2

datax97<-read.csv("x97.csv",header=T)  #测试版 ##2

colnames(datax97)[1] <- 'gene_id' ##1

write.csv(datax97,file="+A56Actin.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################98##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data98<-read.csv('A56Actin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x98<-Reduce(intersect,list(data1$gene_id,
                           data98$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A98<-data1[x98,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A98,file="x98.csv") ##2

datax98<-read.csv("x98.csv",header=T)  #测试版 ##2

colnames(datax98)[1] <- 'gene_id' ##1

write.csv(datax98,file="+A56Actin.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################99##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data99<-read.csv('A57Tropomyosin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x99<-Reduce(intersect,list(data1$gene_id,
                           data99$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A99<-data1[x99,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A99,file="x99.csv") ##2

datax99<-read.csv("x99.csv",header=T)  #测试版 ##2

colnames(datax99)[1] <- 'gene_id' ##1

write.csv(datax99,file="+A57Tropomyosin.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################100##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data100<-read.csv('A57Tropomyosin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x100<-Reduce(intersect,list(data1$gene_id,
                           data100$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A100<-data1[x100,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A100,file="x100.csv") ##2

datax100<-read.csv("x100.csv",header=T)  #测试版 ##2

colnames(datax100)[1] <- 'gene_id' ##1

write.csv(datax100,file="+A57Tropomyosin.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################101##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data101<-read.csv('A58Muscle mitochondrial protein (myosin light chain).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x101<-Reduce(intersect,list(data1$gene_id,
                            data101$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A101<-data1[x101,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A101,file="x101.csv") ##2

datax101<-read.csv("x101.csv",header=T)  #测试版 ##2

colnames(datax101)[1] <- 'gene_id' ##1

write.csv(datax101,file="+A58Muscle mitochondrial protein (myosin light chain).B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################102##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data102<-read.csv('B1TNF-related apoptosis-inducing ligand.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x102<-Reduce(intersect,list(data1$gene_id,
                            data102$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A102<-data1[x102,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A102,file="x102.csv") ##2

datax102<-read.csv("x102.csv",header=T)  #测试版 ##2

colnames(datax102)[1] <- 'gene_id' ##1

write.csv(datax102,file="+B1TNF-related apoptosis-inducing ligand.B1all.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################103##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data103<-read.csv('B1TNF-related apoptosis-inducing ligand.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x103<-Reduce(intersect,list(data1$gene_id,
                            data103$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A103<-data1[x103,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A103,file="x103.csv") ##2

datax103<-read.csv("x103.csv",header=T)  #测试版 ##2

colnames(datax103)[1] <- 'gene_id' ##1

write.csv(datax103,file="+B1TNF-related apoptosis-inducing ligand.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################104##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data104<-read.csv('B2Fas ligand (FasL).B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x104<-Reduce(intersect,list(data1$gene_id,
                            data104$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A104<-data1[x104,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A104,file="x104.csv") ##2

datax104<-read.csv("x104.csv",header=T)  #测试版 ##2

colnames(datax104)[1] <- 'gene_id' ##1

write.csv(datax104,file="+B2Fas ligand (FasL).B1all.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################105##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data105<-read.csv('B2Fas ligand (FasL).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x105<-Reduce(intersect,list(data1$gene_id,
                            data105$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A105<-data1[x105,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A105,file="x105.csv") ##2

datax105<-read.csv("x105.csv",header=T)  #测试版 ##2

colnames(datax105)[1] <- 'gene_id' ##1

write.csv(datax105,file="+B2Fas ligand (FasL).B1someright.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################106##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data106<-read.csv('B3Tumor Necrosis Factor-alpha.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x106<-Reduce(intersect,list(data1$gene_id,
                            data106$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A106<-data1[x106,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A106,file="x106.csv") ##2

datax106<-read.csv("x106.csv",header=T)  #测试版 ##2

colnames(datax106)[1] <- 'gene_id' ##1

write.csv(datax106,file="+B3Tumor Necrosis Factor-alpha.B1all.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################107##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data107<-read.csv('B3Tumor Necrosis Factor-alpha.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x107<-Reduce(intersect,list(data1$gene_id,
                            data107$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A107<-data1[x107,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A107,file="x107.csv") ##2

datax107<-read.csv("x107.csv",header=T)  #测试版 ##2

colnames(datax107)[1] <- 'gene_id' ##1

write.csv(datax107,file="+B3Tumor Necrosis Factor-alpha.B1someright.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################108##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data108<-read.csv('B3Tumor Necrosis Factor-alpha.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x108<-Reduce(intersect,list(data1$gene_id,
                            data108$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A108<-data1[x108,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A108,file="x108.csv") ##2

datax108<-read.csv("x108.csv",header=T)  #测试版 ##2

colnames(datax108)[1] <- 'gene_id' ##1

write.csv(datax108,file="+B3Tumor Necrosis Factor-alpha.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################109##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data109<-read.csv('B4Interleukin-1.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x109<-Reduce(intersect,list(data1$gene_id,
                            data109$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A109<-data1[x109,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A109,file="x109.csv") ##2

datax109<-read.csv("x109.csv",header=T)  #测试版 ##2

colnames(datax109)[1] <- 'gene_id' ##1

write.csv(datax109,file="+B4Interleukin-1.B1all.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################110##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data110<-read.csv('B4Interleukin-1.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x110<-Reduce(intersect,list(data1$gene_id,
                            data110$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A110<-data1[x110,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A110,file="x110.csv") ##2

datax110<-read.csv("x110.csv",header=T)  #测试版 ##2

colnames(datax110)[1] <- 'gene_id' ##1

write.csv(datax110,file="+B4Interleukin-1.B1someright.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################111##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data111<-read.csv('B5Interferon gamma.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x111<-Reduce(intersect,list(data1$gene_id,
                            data111$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A111<-data1[x111,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A111,file="x111.csv") ##2

datax111<-read.csv("x111.csv",header=T)  #测试版 ##2

colnames(datax111)[1] <- 'gene_id' ##1

write.csv(datax111,file="+B5Interferon gamma.B1all.csv",row.names = FALSE) ## 2

###




###数据导出(命名)csv###################################################112##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data112<-read.csv('B5Interferon gamma.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x112<-Reduce(intersect,list(data1$gene_id,
                            data112$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A112<-data1[x112,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A112,file="x112.csv") ##2

datax112<-read.csv("x112.csv",header=T)  #测试版 ##2

colnames(datax112)[1] <- 'gene_id' ##1

write.csv(datax112,file="+B5Interferon gamma.B1someright.csv",row.names = FALSE) ## 2

###




###数据导出(命名)csv###################################################113##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data113<-read.csv('B6TWEAK.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x113<-Reduce(intersect,list(data1$gene_id,
                            data113$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A113<-data1[x113,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A113,file="x113.csv") ##2

datax113<-read.csv("x113.csv",header=T)  #测试版 ##2

colnames(datax113)[1] <- 'gene_id' ##1

write.csv(datax113,file="+B6TWEAK.B1all.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################114##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data114<-read.csv('B6TWEAK.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x114<-Reduce(intersect,list(data1$gene_id,
                            data114$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A114<-data1[x114,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A114,file="x114.csv") ##2

datax114<-read.csv("x114.csv",header=T)  #测试版 ##2

colnames(datax114)[1] <- 'gene_id' ##1

write.csv(datax114,file="+B6TWEAK.B1someright.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################115##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data115<-read.csv('B7Perforin.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x115<-Reduce(intersect,list(data1$gene_id,
                            data115$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A115<-data1[x115,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A115,file="x115.csv") ##2

datax115<-read.csv("x115.csv",header=T)  #测试版 ##2

colnames(datax115)[1] <- 'gene_id' ##1

write.csv(datax115,file="+B7Perforin.B1all.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################116##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data116<-read.csv('B7Perforin.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x116<-Reduce(intersect,list(data1$gene_id,
                            data116$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A116<-data1[x116,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A116,file="x116.csv") ##2

datax116<-read.csv("x116.csv",header=T)  #测试版 ##2

colnames(datax116)[1] <- 'gene_id' ##1

write.csv(datax116,file="+B7Perforin.B1someright.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################117##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data117<-read.csv('B8Granzyme.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x117<-Reduce(intersect,list(data1$gene_id,
                            data117$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A117<-data1[x117,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A117,file="x117.csv") ##2

datax117<-read.csv("x117.csv",header=T)  #测试版 ##2

colnames(datax117)[1] <- 'gene_id' ##1

write.csv(datax117,file="+B8Granzyme.B1all.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################117##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data117<-read.csv('B8Granzyme.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x117<-Reduce(intersect,list(data1$gene_id,
                            data117$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A117<-data1[x117,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A117,file="x117.csv") ##2

datax117<-read.csv("x117.csv",header=T)  #测试版 ##2

colnames(datax117)[1] <- 'gene_id' ##1

write.csv(datax117,file="+B8Granzyme.B1someright.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################118##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data118<-read.csv('B9MICA.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x118<-Reduce(intersect,list(data1$gene_id,
                            data118$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A118<-data1[x118,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A118,file="x118.csv") ##2

datax118<-read.csv("x118.csv",header=T)  #测试版 ##2

colnames(datax118)[1] <- 'gene_id' ##1

write.csv(datax118,file="+B9MICA.B1all.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################119##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data119<-read.csv('B9MICA.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x119<-Reduce(intersect,list(data1$gene_id,
                            data119$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A119<-data1[x119,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A119,file="x119.csv") ##2

datax119<-read.csv("x119.csv",header=T)  #测试版 ##2

colnames(datax119)[1] <- 'gene_id' ##1

write.csv(datax119,file="+B9MICA.B1someright.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################120##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data120<-read.csv('B10MICB.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x120<-Reduce(intersect,list(data1$gene_id,
                            data120$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A120<-data1[x120,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A120,file="x120.csv") ##2

datax120<-read.csv("x120.csv",header=T)  #测试版 ##2

colnames(datax120)[1] <- 'gene_id' ##1

write.csv(datax120,file="+B10MICB.B1all.csv",row.names = FALSE) ## 2

###

###数据导出(命名)csv###################################################121##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data121<-read.csv('B10MICB.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x121<-Reduce(intersect,list(data1$gene_id,
                            data121$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A121<-data1[x121,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A121,file="x121.csv") ##2

datax121<-read.csv("x121.csv",header=T)  #测试版 ##2

colnames(datax121)[1] <- 'gene_id' ##1

write.csv(datax121,file="+B10MICB.B1someright.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################122##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data122<-read.csv('B11immunoglobulins.Ig.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x122<-Reduce(intersect,list(data1$gene_id,
                            data122$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A122<-data1[x122,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A122,file="x122.csv") ##2

datax122<-read.csv("x122.csv",header=T)  #测试版 ##2

colnames(datax122)[1] <- 'gene_id' ##1

write.csv(datax122,file="+B11immunoglobulins.Ig.B1all.csv",row.names = FALSE) ## 2

###



###数据导出(命名)csv###################################################123##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data123<-read.csv('B11immunoglobulins.Ig.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x123<-Reduce(intersect,list(data1$gene_id,
                            data123$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A123<-data1[x123,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A123,file="x123.csv") ##2

datax123<-read.csv("x123.csv",header=T)  #测试版 ##2

colnames(datax123)[1] <- 'gene_id' ##1

write.csv(datax123,file="+B11immunoglobulins.Ig.B1someright.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################124##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data124<-read.csv('B12lysozyme.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x124<-Reduce(intersect,list(data1$gene_id,
                            data124$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A124<-data1[x124,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A124,file="x124.csv") ##2

datax124<-read.csv("x124.csv",header=T)  #测试版 ##2

colnames(datax124)[1] <- 'gene_id' ##1

write.csv(datax124,file="+B12lysozyme.B1all.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################125##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data125<-read.csv('B12lysozyme.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x125<-Reduce(intersect,list(data1$gene_id,
                            data125$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A125<-data1[x125,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A125,file="x125.csv") ##2

datax125<-read.csv("x125.csv",header=T)  #测试版 ##2

colnames(datax125)[1] <- 'gene_id' ##1

write.csv(datax125,file="+B12lysozyme.B1someright.csv",row.names = FALSE) ## 2

###


###数据导出(命名)csv###################################################126##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data126<-read.csv('B13complement.B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x126<-Reduce(intersect,list(data1$gene_id,
                            data126$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A126<-data1[x126,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###

write.csv(A126,file="x126.csv") ##2

datax126<-read.csv("x126.csv",header=T)  #测试版 ##2

colnames(datax126)[1] <- 'gene_id' ##1

write.csv(datax126,file="+B13complement.B1all.csv",row.names = FALSE) ## 2

###126


###数据导出(命名)csv###################################################127##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data127<-read.csv('B13complement.B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x127<-Reduce(intersect,list(data1$gene_id,
                            data127$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A127<-data1[x127,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###127

write.csv(A127,file="x127.csv") ##2

datax127<-read.csv("x127.csv",header=T)  #测试版 ##2

colnames(datax127)[1] <- 'gene_id' ##1

write.csv(datax127,file="+B13complement.B1someright.csv",row.names = FALSE) ## 2

###127



###数据导出(命名)csv###################################################128##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data128<-read.csv('B14Interleukin-8 (IL-8)B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x128<-Reduce(intersect,list(data1$gene_id,
                            data128$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A128<-data1[x128,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###128

write.csv(A128,file="x128.csv") ##2

datax128<-read.csv("x128.csv",header=T)  #测试版 ##2

colnames(datax128)[1] <- 'gene_id' ##1

write.csv(datax128,file="+B14Interleukin-8 (IL-8)B1all.csv",row.names = FALSE) ## 2

###128



###数据导出(命名)csv###################################################129##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data129<-read.csv('B14Interleukin-8 (IL-8)B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x129<-Reduce(intersect,list(data1$gene_id,
                            data129$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A129<-data1[x129,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###129

write.csv(A129,file="x129.csv") ##2

datax129<-read.csv("x129.csv",header=T)  #测试版 ##2

colnames(datax129)[1] <- 'gene_id' ##1

write.csv(datax129,file="+B14Interleukin-8 (IL-8)B1someright.csv",row.names = FALSE) ## 2

###129



###数据导出(命名)csv###################################################130##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data130<-read.csv('B15TNF-related apoptosis-inducing ligand (TRAIL).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x130<-Reduce(intersect,list(data1$gene_id,
                            data130$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A130<-data1[x130,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###130

write.csv(A130,file="x130.csv") ##2

datax130<-read.csv("x130.csv",header=T)  #测试版 ##2

colnames(datax130)[1] <- 'gene_id' ##1

write.csv(datax130,file="+B15TNF-related apoptosis-inducing ligand (TRAIL).B1someright.csv",row.names = FALSE) ## 2

###130



###数据导出(命名)csv###################################################131##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data131<-read.csv('B16Interleukin-1 (IL-1).B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x131<-Reduce(intersect,list(data1$gene_id,
                            data131$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A131<-data1[x131,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###131

write.csv(A131,file="x131.csv") ##2

datax131<-read.csv("x131.csv",header=T)  #测试版 ##2

colnames(datax131)[1] <- 'gene_id' ##1

write.csv(datax131,file="+B16Interleukin-1 (IL-1).B1all.csv",row.names = FALSE) ## 2

###131




###数据导出(命名)csv###################################################132##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data132<-read.csv('B16Interleukin-1 (IL-1).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x132<-Reduce(intersect,list(data1$gene_id,
                            data132$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A132<-data1[x132,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###132

write.csv(A132,file="x132.csv") ##2

datax132<-read.csv("x132.csv",header=T)  #测试版 ##2

colnames(datax132)[1] <- 'gene_id' ##1

write.csv(datax132,file="+B16Interleukin-1 (IL-1).B1someright.csv",row.names = FALSE) ## 2

###132


###数据导出(命名)csv###################################################133##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data133<-read.csv('B17Tumor necrosis factor (TNF).B1all.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x133<-Reduce(intersect,list(data1$gene_id,
                            data133$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A133<-data1[x133,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###133

write.csv(A133,file="x133.csv") ##2

datax133<-read.csv("x133.csv",header=T)  #测试版 ##2

colnames(datax133)[1] <- 'gene_id' ##1

write.csv(datax133,file="+B17Tumor necrosis factor (TNF).B1all.csv",row.names = FALSE) ## 2

###133


###数据导出(命名)csv###################################################134##$
#########
data1<-read.csv("DEG.csv",header=T)  #测试版

colnames(data1)[1] <- 'gene_id'

data134<-read.csv('B17Tumor necrosis factor (TNF).B1someright.csv',header=T)##$$DATA和文档名字要改SSSSSS 1 MMMMMMM 1

x134<-Reduce(intersect,list(data1$gene_id,
                            data134$gene_id))########类型的共同项目##$$两个X后需要要改SSSSSSS 2

rownames(data1)<-data1$gene_id

data1[c(1,1000000000),]

A134<-data1[x134,-1]##$$文档名字要改 X后序号要改SSSSSSS 2

###134

write.csv(A134,file="x134.csv") ##2

datax134<-read.csv("x134.csv",header=T)  #测试版 ##2

colnames(datax134)[1] <- 'gene_id' ##1

write.csv(datax134,file="+B17Tumor necrosis factor (TNF).B1someright.csv",row.names = FALSE) ## 2

###134
#############################ZHEN完结撒花#####################史诗