
### Get the parameters
parser = argparse::ArgumentParser(description = "Script for sankey plot")
parser$add_argument('-i', '--input', help = 'input file name')
opts = parser$parse_args()

library(tidyverse)
library(navdata)
library(networkD3)
library(webshot)
library(data.table)
library(stringi)

# 读取数据，指定
file <- opts$input
data <- read.table(file, sep='\t', header=TRUE)
colnames(data) <- c('source', 'target', 'Value')

# 使用 data.table 对数据进行排序
data <- data.table(data)
data <- data.table(data[order(Value, decreasing=T)])

# 创建节点对应的 id 
node <- unique(c(data$target, data$source))
id <- rep(1:length(node))
node <- data.table(id=id-1, node=node)
for (i in data$source){data$source[data$source==i] <- node$id[node$node==i]}
for (i in data$target){data$target[data$target==i] <- node$id[node$node==i]}
data$source <- as.integer(data$source)
data$target <- as.integer(data$target)

# 对 node 分组，指定 node 颜色
node[stri_detect_regex(node, '^sc'), `:=`(group = 'sc', color='#00ced1')]
node[stri_detect_regex(node, '^Da'), `:=`(group = 'Da', color='#b5cde9')]

colors <- paste(unique(node[, .(group, color)])$color, collapse = '", "')
colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')

p=sankeyNetwork(
  Links=data, Nodes=node, 
  Source="source", Target="target", NodeID="node", Value="Value", 
  NodeGroup="group", colourScale=colorJS,
  nodeWidth=30, nodePadding=58, fontSize=20, sinksRight=F,
  height=800
)
p
saveNetwork(p, paste(file, '.html'))
webshot(paste(file, '.html') , paste(file, ".pdf"))
