# Librairies ####
library(tidyverse)
library(gplots)

# Functions ####


# Parameters ####

# Main ####


# Create file with number of over and under expressed per conditions ----
all_conds = list.dirs("DEG")[-1]

for (i in all_conds){
  actual_cond = strsplit(i, "/")[[1]][3]
  DEG_table = read_tsv(file = paste0(i, "/DGE_result_",actual_cond,"_FC_and_FDR_filter.tsv"),
                       col_names = T)
  nbr_overexpressed = length(which(DEG_table$logFC > 0))
  nbr_underexpressed = length(which(DEG_table$logFC < 0))
  total = nbr_overexpressed + nbr_underexpressed
  df = data.frame(nbr_of_overexpressed = nbr_overexpressed,
                  nbr_of_underexpressed = nbr_underexpressed,
                  total = total)
  write.table(df, file = paste0(i,"/nbr_of_under_and_over_expressed.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
} 

# create matrix for heatmap ----

cond_names = c("0dpa","1dpa","2dpa","3dpa","5dpa", "NA")

cmat = matrix(0, nrow = 6, ncol = 6)

rownames(cmat) = cond_names
colnames(cmat) = cond_names


cmat


# Create overexpressed
for (x in 1:5) {
  for (y in (x+1):6) {
    cond = paste0(cond_names[x],"vs",cond_names[y])
    nbr_under_over = read_tsv(file = paste0("DEG/",cond, "/nbr_of_under_and_over_expressed.txt"),
                              col_names = TRUE)
    cmat[x,y] = nbr_under_over$nbr_of_overexpressed
    cmat[y,x] = nbr_under_over$nbr_of_underexpressed
  }
}

# Create underexpressed
for (x in 1:5) {
  for (y in (x+1):6) {
    cond = paste0(cond_names[x],"vs",cond_names[y])
    nbr_under_over = read_tsv(file = paste0("DEG/",cond, "/nbr_of_under_and_over_expressed.txt"),
                              col_names = TRUE)
    cmat[y,x] = nbr_under_over$nbr_of_underexpressed
    cmat[x,y] = nbr_under_over$nbr_of_underexpressed
  }
}


colors = c(seq(0,1500,length=10),seq(1501,3000,length=10), seq(3000,9000,length=10))
#mypalette = colorRampPalette(c("grey","#1fe1f2","#1f97f2"))(n = 29)

mypalette = colorRampPalette(c("grey","#f27a1f","#f21f1f"))(n = 29)
mypalette = colorRampPalette(c("grey","lightblue1","dodgerblue"))(n = 29)

heatmap.2(cmat, cellnote = cmat,notecol = "black", col=mypalette, key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=1, cexRow = 1,
          main = "Number of Up and Down genes", margins = c(8,8),
          labRow = rownames(cmat), labCol = colnames(cmat),
          dendrogram = 'none', symkey = F, symm = F, symbreaks = F,
          Rowv = F, Colv = F,
          scale = "none", xlab = "", ylab = "", offsetRow=-30)


conds1 = c(rep(rownames(cmat)[1], 5),
           rep(rownames(cmat)[2], 4),
           rep(rownames(cmat)[3], 3),
           rep(rownames(cmat)[4], 2),
           rep(rownames(cmat)[5], 1))

conds2 = c(rownames(cmat)[2:6],
           rownames(cmat)[3:6],
           rownames(cmat)[4:6],
           rownames(cmat)[5:6],
           rownames(cmat)[6])


val2 = c(cmat[2:6,1], cmat[3:6,2], cmat[4:6,3], cmat[5:6,4], cmat[6,5])
val1 = c(cmat[1,2:6], cmat[2,3:6], cmat[3,4:6], cmat[4,5:6], cmat[5,6])

tib2 = tibble(
  g1 = conds1,
  g2 = conds2,
  val1 = val1,
  val2 = val2
)

tib2_a = asymmetrise(tib2, g1, g2)
tib2_a[31:36,3:4] = 0

ggplot(tib2_a, aes(x = g1, y = g2)) +
  geom_asymmat(aes(fill_tl = val1, fill_br = val2, fill_diag = rep(0,36))) +
  scale_fill_tl_gradient(low = "lightpink", high = "tomato") +
  scale_fill_br_gradient(low = "lightblue1", high = "dodgerblue") +
  scale_fill_diag_gradient(low = "white", high = "grey",
                           guide = guide_none())+
  labs(
    fill_tl = "Up-regulated genes",
    fill_br = "Down-regulated genes"
  )+
  theme_minimal()+
  theme(axis.title = element_blank())