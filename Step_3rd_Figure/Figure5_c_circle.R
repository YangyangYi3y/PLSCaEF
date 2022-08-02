library(circlize)
library(reshape2)
dev.off()
dev.new()
datapath = 'D:\\projects\\ABCD\\EF\\PLS\\results\\ef_task23\\combat\\2fold\\comp1\\mat.txt'
mat <- read.table(datapath,sep='')
mat=apply(mat,2,as.numeric)
mat[lower.tri(mat)] <- 0
rn <- c('Aud','CO','CP','DMN','DAN','FP','None','RSP','SM','SMla','Sal','VAN','Vis','Sub')
rownames(mat) <- rn
colnames(mat) <- rn
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df
grid.col = c(Aud = "#FF0066", CO = "#6e2998", CP = "#CC9966", DMN = "#FF33FF", DAN = "#FF0033", 
             FP = "#98E98A", None = "grey", RSP= "#FF9999", SM = "#FF3366",
             SMlat = "#FF9933", Sal = "#996666", VAN= "#339999", Vis = "#33CC99", Sub = "#333333")
circos.par(gap.after = c(rep(3, nrow(mat))),start.degree = 90)
#col_mat = rand_color(length(mat), transparency = 0)
col_mat[mat < 0.02] = "#00000000"
dim(col_mat) = dim(mat)  # to make sure it is a matrix
chordDiagram(mat, order = rn,col = "red",link.visible = mat >0.1,link.sort = TRUE,grid.col = grid.col,link.decreasing = TRUE, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.06))
circos.clear()
