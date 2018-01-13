#!/usr/bin/Rscript
library(getopt)
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
# 一般就是4列，第一列为字符串，第二列为简写，第三列值分别为0（无参数后面可以不跟参数）、1（后面需要跟参数）、2（可选可不选），第四列为数据类型
# character logical integer double
spec = matrix(c(
  'input_file', 'i', 1, "character",
  'help'  , 'h', 0, "logical",
  'count' , 'c', 1, "character", # 1:17,2:28,3:28,4:28 分别指出每组分类（1,2,3,4）含有多少个sample
  'plot_name' , 'n', 1, "character", # you must point the name of the plot
  'output_file' , 'o' , 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
# if ( is.null(opt$output_file) ) { opt$output_file = opt$input_file }

library(ggplot2)
library(reshape)
#增加分组信息(genus)
data<- read.table(paste(opt$input_file),sep='\t',header=T,check.names=F)
rank <-data$tax_name
row=nrow(data)
b<-strsplit(paste(opt$count),split = ",")
i=1
group=c()
for (tmp in b){
  #print(tmp)
  t<-strsplit(tmp,split = ":")
  #print(t)
}
for (tmp in t){
  group<-c(group,c(rep(t[[i]][1],row*as.numeric(t[[i]][2]))))
  #print(a)
  i=i+1
}

#print(group)
#group <-c(rep(1,row*17),rep(2,row*28),rep(3,row*30),rep(4,row*28))
data <-melt.data.frame(data,id="tax_name")
data <- cbind(data,group)
#设置固定的tax_name顺序
data$tax_name = factor(data$tax_name, levels=rank) 
#colours <- c("#FFFFFF","#F4F8FC","#E8F1F9","#DDE9F5","#D1E2F2","#C6DBEF","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b","#08306b")
#设置背景模板
theme_set(theme_bw())
ggplot(data, aes(x = variable,y=value,fill=tax_name) )+ geom_bar(stat = "identity",colour="black",size=0.2) +#geom_histogram(binwidth=10)
  #scale_fill_brewer(palette="Blues")+ #Greens
  theme(axis.text.x = element_text(angle = 90,hjust = .5, vjust = .5,size = 5)) +
  facet_grid(. ~ group,scales = "free_x",space = "free_x")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold",size = rel(1.5)),strip.background = element_rect(fill = "lightblue",colour = "black",size = 1))+
  ggtitle(paste(opt$plot_name))+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Samples")+ylab("Relative Abundance(%)")+
  #scale_fill_manual(values=colours[1:(N+1)])
  scale_fill_hue()+
  # change the legend font
  theme(legend.text = element_text( size = 1.5),legend.key.size=unit(0.2,'cm'))+
  theme(legend.position="bottom")
#  theme(legend.position="none") 
ggsave(paste(opt$output_file),width=8,height=8)

## output the legend only (没有解决)
#dir <- dirname(paste(opt$output_file))
#filename <-paste(opt$plot_name)
#legend <- file.path(dir,"legend.pdf")
#require(grid)
#g <- ggplot_gtable(ggplot_build(p))$grobs
#dev.new()
#pdf(legend, width=6, height=3)
##g
#pushViewport(plotViewport(rep(1, 4)))
#grid.draw(g[[24]])
#dev.off()

