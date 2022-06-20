library(vegan)
some_list <- c(
    "Malaya","SM1.Arctoides-1","SM2.Arctoides-2","A20","XH1.Assamensis","BGI-CE-4.fascicularis","BGI-96346.mulatta","BGI.Mulatta-1","CR-5.Mulatta-2","Tibetan-macaque-NO.3"
)

some_list <-c("DRR219369","DRR219370","DRR219371","SRA023855","SRR1024051","SRR1564766","SRR2981114")

seq_to <- seq(1,79) 

some_list <- paste0('C_rhe_',seq_to)
root <- "F:/Work/Now_work/2021_05_18.ben_/15.new_tests/"
root <- "F:/Work/Now_work/2021_05_18.ben_/15.new_tests/Chinese_all/"


input_species <- read.table(file.path(root,"ind_all_species.txt"), sep = " ") 
some_list <- input_species$V1

out_mat <- matrix(nrow = length(some_list), ncol = length(some_list))

for (i in seq(1,length(some_list))) {
    sp1 = some_list[i]
    for (j in seq(i+1,length(some_list))) {
        sp2 = some_list[j]
        to_process = paste0("fst_",sp1,"_",sp2,".density.out")
        
        file_tpath <- file.path(root, "Ind_all/input/output",to_process)
        data_to <- read.table(file_tpath, header = T)
        
        mod <- lm(Fst ~ containsNinteractgenez * number_of_genes ,data = data_to)
        out_mat[i,j] = coef(mod)["containsNinteractgenez"]

    }

}
colnames(out_mat) <- some_list
rownames(out_mat) <- some_list
View(t(out_mat))

write.csv(x = t(out_mat), file = "coef.csv",na = "" ,quote = F)

dat <- read.csv(file.path(root,"Chinese_all","combined.csv"), header = T)
mod <- lm(fst ~ dist_ratio, data = dat)
summary(mod)

#

aDNA <- read.csv(file.path(root,"Ind_all",'Ind_aDNA.csv'), header = F)
mtDNA <- read.csv(file.path(root,"Ind_all",'Ind_mtDNA.csv'), header = F)
out_mat <- read.csv(file.path(root,"Ind_all",'coef.csv'), header = T)
out_mat <- out_mat[,-1]

ad_matrix <- matrix(nrow = length(some_list), ncol = length(some_list))
colnames(ad_matrix) = some_list
rownames(ad_matrix) = some_list
for (i in some_list) {
    for (j in some_list) {
        ad_matrix[i,j] = aDNA[which(aDNA$V1 == i & aDNA$V2 == j),3]
    }
}

mtd_matrix <- matrix(nrow = length(some_list), ncol = length(some_list))
colnames(mtd_matrix) = some_list
rownames(mtd_matrix) = some_list
for (i in some_list) {
    for (j in some_list) {
        mtd_matrix[i,j] = mtDNA[which(mtDNA$V1 == i & mtDNA$V2 == j),3]
    }
}


coef_matrix <- matrix(data = 0,nrow = length(some_list), ncol= length(some_list))
coef_matrix[lower.tri(coef_matrix, diag=FALSE)] <- out_mat[lower.tri(coef_matrix, diag=FALSE)]
coef_matrix[upper.tri(coef_matrix, diag=FALSE)] <- t(coef_matrix)[upper.tri(t(coef_matrix))]

mantel(coef_matrix,mtd_matrix)

mantel.partial(coef_matrix,mtd_matrix,ad_matrix, permutations = 1000)

#
out_mat <- read.csv(file.path(root, "Aureus_all","coef.csv"), header = T)
out_mat <- as.matrix(out_mat[,2:(length(some_list)+1)])

coef_matrix <- matrix(data = 0,nrow = length(some_list), ncol= length(some_list))
coef_matrix[upper.tri(coef_matrix, diag=FALSE)] <- out_mat[lower.tri(coef_matrix, diag=FALSE)]
coef_matrix[lower.tri(coef_matrix, diag=FALSE)] <- t(coef_matrix)[lower.tri(t(coef_matrix))]

ad_matrix <- read.table(file.path(root, "Aureus_all","aur_a.csv"), header = T, sep='\t')
ad_matrix <- as.matrix(ad_matrix[,2:(length(some_list)+1)])

mtd_matrix <- read.table(file.path(root, "Aureus_all","aur_mt.csv"), header = T, sep='\t')
mtd_matrix <- as.matrix(mtd_matrix[,2:(length(some_list)+1)])

ratio_matrix <- mtd_matrix/ad_matrix

mantel(ratio_matrix,coef_matrix)
mantel(coef_matrix,mtd_matrix)

mantel.partial(coef_matrix,mtd_matrix,ad_matrix,permutations=1000)

#
dat$dist_ratio <- dat$mtdist/dat$adist

ratio_matrix <- matrix(data = 0,nrow = 79, ncol = 79)
ratio_matrix[upper.tri(ratio_matrix, diag=FALSE)] <- dat$dist_ratio
ratio_matrix[lower.tri(ratio_matrix, diag=FALSE)] <- t(ratio_matrix)[lower.tri(t(ratio_matrix))]

adist_matrix <- matrix(data = 0,nrow = 79, ncol = 79)
adist_matrix[upper.tri(adist_matrix, diag=FALSE)] <- dat$adist
adist_matrix[lower.tri(adist_matrix, diag=FALSE)] <- t(adist_matrix)[lower.tri(t(adist_matrix))]

mtdist_matrix <- matrix(data = 0,nrow = 79, ncol = 79)
mtdist_matrix[upper.tri(mtdist_matrix, diag=FALSE)] <- dat$mtdist
mtdist_matrix[lower.tri(mtdist_matrix, diag=FALSE)] <- t(mtdist_matrix)[lower.tri(t(mtdist_matrix))]

coef_matrix <- matrix(data = 0,nrow = 79, ncol = 79)
coef_matrix[upper.tri(coef_matrix, diag=FALSE)] <- dat$Coef_containN
coef_matrix[lower.tri(coef_matrix, diag=FALSE)] <- t(coef_matrix)[lower.tri(t(coef_matrix))]

############

mantel(ratio_matrix,coef_matrix)
mantel(coef_matrix,mtdist_matrix)

mantel.partial(coef_matrix,mtdist_matrix,adist_matrix,permutations=1000)
