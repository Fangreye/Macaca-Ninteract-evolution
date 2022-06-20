library("VennDiagram")
library('eulerr')
library('UpSetR')

setwd("F:/Work/Now_work/2021_05_18.ben_/12.euler_plot")
#### fst
chinese_group <- list(
  "Red : Blue" = c('TARS2','TIMMDC1','ATP5O','MRPS36','MRPL52','NDUFB8','MRPL51','ACP1','MRPL11','MRPL49','SURF1','MRPL10','MRPL54','PET100','MRPL4','NDUFA13','ATP5SL','MRPS34','COX4I1'),
  "Red : Purple" = c('MRPL24','ATP5O','MRPS36','CD14','NDUFA2','HARS2','TYMP','BCS1L','SURF1','MRPL10','POLRMT','MRPL54','NDUFA13','SARS2','MRPS12','ATP5SL','MRPL28','MRPS34','EARS2'),
  "Blue : Purple" = c('DARS2','ATP5O','NDUFB2','MRPS18A','ATP5I','LYRM7','CD14','NDUFA2','NDUFAF1','COX5A','MRPL52','TFAM','NDUFB8','MRPL43','TYMP','ACP1','UQCC3','ATP5L','SURF1','MRPS7','UQCR','MRPL54','NDUFA7','MRPL4','MRPS34','COX4I1
  ')
)

indian_group <- list(
  "Orange : Brown" = c("MRPL3","UQCRQ","NDUFB8","LOC716161","MRPL44","MRPL30","MRPL49","FOXRED1","MRPL54","SARS2","MRPS12"),
  'Orange : Red' = c("NDUFS5","LARS2","ATP5MF","TMEM70","CYC1","NDUFB8","SARS2","MRPS12","MRPS34","COX4I1"),
  'Brown : Red' = c("LARS2","ACAD9","ATP5MF","NDUFB8","LOC716161","SARS2","MRPS12","MRPL28","MRPS34","COX4I1")
)

# arctoides <- list(
#   "arctoides : assamensis" = c('DARS2','NDUFAF3','ATP5J2','ATP5I','TMEM70','MRPS16','MRPL40','TYMP','ATP5B','NDUFB3','UQCC3','TTC19','POLRMT','NDUFA13','MRPL28'),     
#   "arctoides : fascicularis" = c('TARS2','MRPL55','NDUFAF3','ATP5J','AARS2','RARS2','CD14','NDUFA2','HARS2','CYC1','MRPS16','TYMP','MRPL51','ATP5B','UQCC3','NDUFA13'),
#   "arctoides : mulata" = c('MRPL55','ATP5J2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','MRPL51','ATP5B','UQCC3','TMEM126B','MRPL34'),
#   "arctoides : thibetana" = c('ATP5I','HARS2','MRPL51','UQCC3','MRPS7','NDUFA7','NDUFA13','MRPL28','MRPS34'),
#   "assamensis : fascicularis" = c('MRPL55','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','MRPL10','MRPL57','POLRMT'),
#   "assamensis : mulata" = c('ATP5F1','MRPL55','DARS2','NDUFAF3','MRPL2','MRPS36','UQCRQ','HARS2','NDUFB1','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','TMEM126B','MRPL45'),     
#   "assamensis : thibetana" = c('MRPL2','ATP5I','MRPS18C','NDUFAF2','HARS2','TMEM70','MRPS16','MRPL43','MRPL40','UQCC3','MRPS2','CCDC56','MRPS7','MRPL38','MRPL57','POLRMT','MRPL20','MRPS34'),
#   "mulata : fascicularis" = c('MRPL37','NDUFAF3','ATP5J2','MRPL2','MRPL1','NDUFC1','MRPS36','LYRM7','CD14','NDUFA2','MRPL22','MRPS16','NDUFB8','MRPL43','USMG5','MRPL51','ATP5B','BCS1L','COX7C','TMEM126B','MRPL10','NDUFA13'),
#   "thibetana : fascicularis" = c('MRPL37','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','CYC1','MRPL43','USMG5','ATP5B','MRPL49','UQCC3','MRPL10','MRPS7','POLRMT','NDUFA13'),       
#   "thibetana : mulata" = c('ATP5F1','MRPL9','MRPL55','MRPL2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','NDUFV1','MRPL49','UQCC3','TMEM126B','MRPL34','NDUFA13')
# )

arctoides <- list(
  "arc : ass" = c('DARS2','NDUFAF3','ATP5J2','ATP5I','TMEM70','MRPS16','MRPL40','TYMP','ATP5B','NDUFB3','UQCC3','TTC19','POLRMT','NDUFA13','MRPL28'),
  "arc : fas" = c('TARS2','MRPL55','NDUFAF3','ATP5J','AARS2','RARS2','CD14','NDUFA2','HARS2','CYC1','MRPS16','TYMP','MRPL51','ATP5B','UQCC3','NDUFA13'),
  "arc : mul" = c('MRPL55','ATP5J2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','MRPL51','ATP5B','UQCC3','TMEM126B','MRPL34'),
  "arc : thi" = c('ATP5I','HARS2','MRPL51','UQCC3','MRPS7','NDUFA7','NDUFA13','MRPL28','MRPS34'),
  "ass : fas" = c('MRPL55','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','MRPL10','MRPL57','POLRMT'),
  "ass : mul" = c('ATP5F1','MRPL55','DARS2','NDUFAF3','MRPL2','MRPS36','UQCRQ','HARS2','NDUFB1','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','TMEM126B','MRPL45'),
  "ass : thi" = c('MRPL2','ATP5I','MRPS18C','NDUFAF2','HARS2','TMEM70','MRPS16','MRPL43','MRPL40','UQCC3','MRPS2','CCDC56','MRPS7','MRPL38','MRPL57','POLRMT','MRPL20','MRPS34'),
  "mul : fas" = c('MRPL37','NDUFAF3','ATP5J2','MRPL2','MRPL1','NDUFC1','MRPS36','LYRM7','CD14','NDUFA2','MRPL22','MRPS16','NDUFB8','MRPL43','USMG5','MRPL51','ATP5B','BCS1L','COX7C','TMEM126B','MRPL10','NDUFA13'),
  "thi : fas" = c('MRPL37','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','CYC1','MRPL43','USMG5','ATP5B','MRPL49','UQCC3','MRPL10','MRPS7','POLRMT','NDUFA13'),
  "thi : mul" = c('ATP5F1','MRPL9','MRPL55','MRPL2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','NDUFV1','MRPL49','UQCC3','TMEM126B','MRPL34','NDUFA13')
)

aureus <- list(
  'aur : thi' = c('MRPS21','NDUFAF3','ATP5MF','MRPL14','UQCRQ','MRPL52','BCS1L','UQCC3','MRPL10','CCDC56'),
  'aur : fas' = c('MRPS21','HIGD1A','NDUFAF3','RARS2','UQCRQ','COX16','NDUFB8','MARS2','NDUFB3','BCS1L','MRPL53','MRPL23','NARS2','TMEM126B','TTC19','ATPAF2','MRPL10','CCDC56','TACO1'),
  'aur : ass' = c('MRPS21','NDUFAF3','ATP5MF','TFB1M','RARS2','MRPS18A','UQCRQ','BCS1L','UQCC3','TTC19','MRPL10','MRPL12','CARS2','NDUFA13'),
  'ass : thi' = c('MRPL20','NDUFAF3','COX19','ATP5MF','RARS2','ATP5ME','MRPL52','CYC1','MRPL21','CARS2','UQCR','NDUFA7','NDUFA13','MRPL28','MRPS34'),
  'fas : thi' = c('NDUFAF3','MRPL2','UQCRQ','NDUFA2','CD14','HARS2','CYC1','BCS1L','COX7C','MRPL49','MRPL10','CCDC56','NDUFA13'),
  'fas : ass' = c('PDC','NDUFAF3','RARS2','MRPL2','UQCRQ','NDUFA2','CD14','HARS2','COX6C','CYC1','BCS1L','COX7C','MRPL49','UQCC3','NARS2','MRPL10','MRPL45','TACO1','NDUFA13')
)

# aureus <- list(
#   "aurea : fascicularis" <- c('MRPS21','HIGD1A','NDUFAF3','RARS2','UQCRQ','COX16','NDUFB8','MARS2','NDUFB3','BCS1L','MRPL53','MRPL23','NARS2','TMEM126B','TTC19','ATPAF2','MRPL10','CCDC56'),
#   "fascicularis : thibetana" <- c('NDUFAF3','MRPL2','UQCRQ','NDUFA2','CD14','HARS2','CYC1','BCS1L','COX7C','MRPL49','UQCC3','NARS2','MRPL10','MRPL45','CCDC56','TACO1','NDUFA13'),
#   "aurea : thibetana" <- c('MRPS21','ATP5PB','NDUFAF3','ATP5MF','MRPL14','UQCRQ','MRPS26','BCS1L','UQCC3','TMEM126B','MRPL10','CCDC56','MRPL12','CARS2','NDUFA13')
# )

chinese_overlap <- get.venn.partitions(chinese_group)
aureus_overlap <- get.venn.partitions(aureus)
indian_overlap  <- get.venn.partitions(indian_group)
arctoides_overlap <- get.venn.partitions(arctoides)

write.csv(View(chinese_overlap))

p <- calculate.overlap(arctoides)

# a <- get.venn.partitions(arctoides)
# a <- a[a$..values.. != character(0),]
# grid.draw(venn.diagram(arctoides, NULL))

te <- fromList(arctoides)
# te <- te[rowSums(te) >= (length(te[4,]) -2),]
pdf("pi_arctoides.pdf", w=8, h=6.0, version="1.4", bg="transparent")

upset(te, sets= rev(names(arctoides)), #nsets = 10, 
      nintersects = NA, 
      text.scale= c(2.5, 2.5, 2, 2, 2.5, 2.5), 
      point.size = 3,
      # queries = list(list(query = intersects, params = list("arc : ass",
      #                                                       "arc : fas"), color = "red", active = T)),
      keep.order = TRUE,
      order.by = 'degree'
) 
dev.off()

# w=8, h=6.0
# w=15, h=8.0
c("thi/ass","fas","aur")
c(2.5, 2.5, 2, 2, 2.5, 2.5)
c(2.5, 2.5, 2, 2, 2.3, 2.5)
text.scale=c(2.5, 2.5, 2, 2, 2.5, 2.5)
c(2.5, 2, 2, 2, 2.5, 2.5)

upset(te)
intersects(te)

aaa(te)

te <- fromList(aureus)
te <- unique(te)
result_list <- rep("",nrow(te))
common_count <- rep(0, nrow(te))
for (i in seq(1, nrow(te))) {
    this <- which(te[i,] == 1)
    text <- Reduce(intersect, aureus[this])
    result_list[i] <- paste(text, collapse = ';')
    common_count[i] <- length(text)
}

write.csv(data.frame(te, 'common_count'= common_count, 'intersect'= result_list), 
          file = "F:/Work/Now_work/2021_05_18.ben_/12.euler_plot/fst_aureus.csv",
          quote = F, row.names = F)



####

z = c()
for (i in seq(1:length(which(aureus_overlap$..count.. > 0)))) {
  z[i] = paste0('"', toString(unlist(aureus_overlap$..values..[which(aureus_overlap$..count.. > 0)][i])), '"' )
}

cols <- sapply(aureus_overlap, is.logical)
num <- length(which(cols))
aureus_overlap[,cols] <- lapply(aureus_overlap[,cols], as.numeric)

a <- data.frame(
  aureus_overlap[1:num][which(aureus_overlap$..count.. > 0),],
  genes = z,
  count = aureus_overlap$..count..[aureus_overlap$..count.. > 0]
)
write.csv(a,"pi_aureus.csv",quote = F, row.names = F)



#### pi
chinese_group <- list(
  'Red' = c('MRPL9','MRPL55','ATP5J2','UQCRQ','ATP5S','NDUFB8','MRPL43','MRPL53','MRPL49','TMEM126B','MRPL10','NDUFA7'),
  'Blue' = c('MRPL55','ATP5J2','HARS2','NDUFB8','MRPL43','MRPL53','MRPL49','TMEM126B','NDUFA7'),
  'Purple' = c('MRPL55','ATP5J2','UQCRQ','HARS2','NDUFB8','MRPL43','MRPL53','NDUFV1','MRPL49','TMEM126B','NDUFA7')
)

indian_group <- list(
  'Orange' = c("MRPS21","MRPL18","MRPL2","UQCRQ","ATP5MPL","NDUFB8","MRPL43","MRPL10"),
  'Brown' = c("MRPS21","MRPL18","MRPL2","UQCRQ","ATP5MPL","MRPL43","NDUFS1","UQCC3","MRPL10","MRPS7"),
  'Red' = c("ATP5MF","MRPL2","UQCRQ","ATP5MPL","TMEM70","NDUFB8","LOC716161","SARS2","MRPS12","MRPL28","MRPS34","EARS2","COX4I1")
)

arctoides <- list(
  'arc' = c('MRPL55','ATP5I','TYMP','MRPL51','MRPS34'),
  'fas' = c('MRPL37','MRPL2','CD14','NDUFA2','HARS2','CYC1','MRPL43','USMG5','MRPL49','UQCC3','MRPL10','POLRMT','NDUFA13'),
  'mul' = c('MRPL37','ATP5F1','MRPL9','MRPL55','ATP5J2','MRPL2','CD14','NDUFA2','HARS2','NDUFB1','NDUFB8','MRPL43','USMG5','NDUFV1','TMEM126B','MRPL50','MRPL34'),   
  'thi' = c('MRPL2','MRPL49','UQCC3','SURF1','CCDC56'),
  'ass' = c('NDUFS5','MRPS21','MRPL55','DARS2','MRPL2','UQCRQ','MRPL15','TMEM70','MRPS16','MRPL43','MRPS26','NDUFB3','BCS1L','UQCC3','TTC19','POLRMT')
)

aureus <- list(
  'aur' = c('MRPS21','NDUFAF3','ATP5MF','UQCRQ','MRPL43','MRPL51','BCS1L','UQCC3','MRPL10','CCDC56','CARS2'),
  'fas' = c('NDUFAF3','MRPL18','RARS2','MRPL2','NDUFAF2','MRPS36','MRPS27','NDUFA2','CD14','NDUFAF6','CYC1','BCS1L','MRPL53','COX7C','MRPL49','NARS2','TMEM126B','MRPL10','CCDC56','NDUFA13','NDUFB10'),
  'thi' = c('MRPL18','MRPL2','NDUFA2','CD14','MRPL49','LOC722212','MRPL10','CCDC56','NDUFB10'),
  'ass' = c('ATP5MF','MRPL36','NDUFS6','UQCC3','MRPL12','CARS2','POLRMT','NDUFA13','MRPS34','NDUFB10')
)

chinese_overlap <- get.venn.partitions(chinese_group)
aureus_overlap <- get.venn.partitions(aureus)
indian_overlap  <- get.venn.partitions(indian_group)
arctoides_overlap <- get.venn.partitions(arctoides)

as.data.frame(arctoides_overlap$..values.., )

write.csv(View(arctoides_overlap))

grid.draw(venn.diagram(a, NULL))

#####################s
whole_group <- list(
    "Red_China : Blue" = c('TARS2','TIMMDC1','ATP5O','MRPS36','MRPL52','NDUFB8','MRPL51','ACP1','MRPL11','MRPL49','SURF1','MRPL10','MRPL54','PET100','MRPL4','NDUFA13','ATP5SL','MRPS34','COX4I1'),
    "Red_China : Purple" = c('MRPL24','ATP5O','MRPS36','CD14','NDUFA2','HARS2','TYMP','BCS1L','SURF1','MRPL10','POLRMT','MRPL54','NDUFA13','SARS2','MRPS12','ATP5SL','MRPL28','MRPS34','EARS2'),
    "Blue : Purple" = c('DARS2','ATP5O','NDUFB2','MRPS18A','ATP5I','LYRM7','CD14','NDUFA2','NDUFAF1','COX5A','MRPL52','TFAM','NDUFB8','MRPL43','TYMP','ACP1','UQCC3','ATP5L','SURF1','MRPS7','UQCR','MRPL54','NDUFA7','MRPL4','MRPS34','COX4I1
  '),
    "Orange : Brown" = c("MRPL3","UQCRQ","NDUFB8","LOC716161","MRPL44","MRPL30","MRPL49","FOXRED1","MRPL54","SARS2","MRPS12"),
    'Orange : Red_India' = c("NDUFS5","LARS2","ATP5MF","TMEM70","CYC1","NDUFB8","SARS2","MRPS12","MRPS34","COX4I1"),
    'Brown : Red_India' = c("LARS2","ACAD9","ATP5MF","NDUFB8","LOC716161","SARS2","MRPS12","MRPL28","MRPS34","COX4I1"),
    "arc : ass" = c('DARS2','NDUFAF3','ATP5J2','ATP5I','TMEM70','MRPS16','MRPL40','TYMP','ATP5B','NDUFB3','UQCC3','TTC19','POLRMT','NDUFA13','MRPL28'),
    "arc : fas" = c('TARS2','MRPL55','NDUFAF3','ATP5J','AARS2','RARS2','CD14','NDUFA2','HARS2','CYC1','MRPS16','TYMP','MRPL51','ATP5B','UQCC3','NDUFA13'),
    "arc : mul" = c('MRPL55','ATP5J2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','MRPL51','ATP5B','UQCC3','TMEM126B','MRPL34'),
    "arc : thi" = c('ATP5I','HARS2','MRPL51','UQCC3','MRPS7','NDUFA7','NDUFA13','MRPL28','MRPS34'),
    "ass : fas" = c('MRPL55','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','MRPL10','MRPL57','POLRMT'),
    "ass : mul" = c('ATP5F1','MRPL55','DARS2','NDUFAF3','MRPL2','MRPS36','UQCRQ','HARS2','NDUFB1','MRPS16','MRPL43','USMG5','ATP5B','UQCC3','TMEM126B','MRPL45'),
    "ass : thi" = c('MRPL2','ATP5I','MRPS18C','NDUFAF2','HARS2','TMEM70','MRPS16','MRPL43','MRPL40','UQCC3','MRPS2','CCDC56','MRPS7','MRPL38','MRPL57','POLRMT','MRPL20','MRPS34'),
    "mul : fas" = c('MRPL37','NDUFAF3','ATP5J2','MRPL2','MRPL1','NDUFC1','MRPS36','LYRM7','CD14','NDUFA2','MRPL22','MRPS16','NDUFB8','MRPL43','USMG5','MRPL51','ATP5B','BCS1L','COX7C','TMEM126B','MRPL10','NDUFA13'),
    "thi : fas" = c('MRPL37','NDUFAF3','MRPL2','CD14','NDUFA2','HARS2','CYC1','MRPL43','USMG5','ATP5B','MRPL49','UQCC3','MRPL10','MRPS7','POLRMT','NDUFA13'),
    "thi : mul" = c('ATP5F1','MRPL9','MRPL55','MRPL2','CD14','NDUFA2','HARS2','NDUFB8','MRPL43','NDUFV1','MRPL49','UQCC3','TMEM126B','MRPL34','NDUFA13'),
    'aur : thi' = c('MRPS21','NDUFAF3','ATP5MF','MRPL14','UQCRQ','MRPL52','BCS1L','UQCC3','MRPL10','CCDC56'),
    'aur : fas' = c('MRPS21','HIGD1A','NDUFAF3','RARS2','UQCRQ','COX16','NDUFB8','MARS2','NDUFB3','BCS1L','MRPL53','MRPL23','NARS2','TMEM126B','TTC19','ATPAF2','MRPL10','CCDC56','TACO1'),
    'aur : ass' = c('MRPS21','NDUFAF3','ATP5MF','TFB1M','RARS2','MRPS18A','UQCRQ','BCS1L','UQCC3','TTC19','MRPL10','MRPL12','CARS2','NDUFA13'),
    'ass : thi' = c('MRPL20','NDUFAF3','COX19','ATP5MF','RARS2','ATP5ME','MRPL52','CYC1','MRPL21','CARS2','UQCR','NDUFA7','NDUFA13','MRPL28','MRPS34'),
    'fas : thi' = c('NDUFAF3','MRPL2','UQCRQ','NDUFA2','CD14','HARS2','CYC1','BCS1L','COX7C','MRPL49','MRPL10','CCDC56','NDUFA13'),
    'fas : ass' = c('PDC','NDUFAF3','RARS2','MRPL2','UQCRQ','NDUFA2','CD14','HARS2','COX6C','CYC1','BCS1L','COX7C','MRPL49','UQCC3','NARS2','MRPL10','MRPL45','TACO1','NDUFA13')
)

te <- fromList(whole_group)
# te <- te[rowSums(te) >= (length(te[4,]) -2),]
pdf("pi_whole.pdf", w=12, h=10.0, version="1.4", bg="transparent")

upset(te, sets= rev(names(whole_group)), #nsets = 10, 
      nintersects = NA, 
      text.scale= c(2.5, 2, 2, 2, 1.6, 2.5), 
      # 1.3
      point.size = 3,
      # queries = list(list(query = intersects, params = list("arc : ass",
      #                                                       "arc : fas"), color = "red", active = T)),
      keep.order = TRUE,
      order.by = 'degree'
) 
dev.off()

a <- overlapGroups(whole_group)

overlap(whole_group[[1]],whole_group[[2]])

te <- fromList(whole_group)
te <- unique(te)
result_list <- rep("",nrow(te))
common_count <- rep(0, nrow(te))
for (i in seq(1, nrow(te))) {
    this <- which(te[i,] == 1)
    text <- Reduce(intersect, whole_group[this])
    result_list[i] <- paste(text, collapse = ';')
    common_count[i] <- length(text)
}

write.csv(data.frame(te, 'common_count'= common_count, 'intersect'= result_list), 
          file = "F:/Work/Now_work/2021_05_18.ben_/12.euler_plot/pi_whole.csv",
          quote = F, row.names = F)


##################
whole_group <- list(
    'Red_China' = c('MRPL9','MRPL55','ATP5J2','UQCRQ','ATP5S','NDUFB8','MRPL43','MRPL53','MRPL49','TMEM126B','MRPL10','NDUFA7'),
    'Blue' = c('MRPL55','ATP5J2','HARS2','NDUFB8','MRPL43','MRPL53','MRPL49','TMEM126B','NDUFA7'),
    'Purple' = c('MRPL55','ATP5J2','UQCRQ','HARS2','NDUFB8','MRPL43','MRPL53','NDUFV1','MRPL49','TMEM126B','NDUFA7'),
    'Orange' = c("MRPS21","MRPL18","MRPL2","UQCRQ","ATP5MPL","NDUFB8","MRPL43","MRPL10"),
    'Brown' = c("MRPS21","MRPL18","MRPL2","UQCRQ","ATP5MPL","MRPL43","NDUFS1","UQCC3","MRPL10","MRPS7"),
    'Red_India' = c("ATP5MF","MRPL2","UQCRQ","ATP5MPL","TMEM70","NDUFB8","LOC716161","SARS2","MRPS12","MRPL28","MRPS34","EARS2","COX4I1"),
    'arc' = c('MRPL55','ATP5I','TYMP','MRPL51','MRPS34'),
    'fas_1' = c('MRPL37','MRPL2','CD14','NDUFA2','HARS2','CYC1','MRPL43','USMG5','MRPL49','UQCC3','MRPL10','POLRMT','NDUFA13'),
    'mul' = c('MRPL37','ATP5F1','MRPL9','MRPL55','ATP5J2','MRPL2','CD14','NDUFA2','HARS2','NDUFB1','NDUFB8','MRPL43','USMG5','NDUFV1','TMEM126B','MRPL50','MRPL34'),   
    'thi_1' = c('MRPL2','MRPL49','UQCC3','SURF1','CCDC56'),
    'ass_1' = c('NDUFS5','MRPS21','MRPL55','DARS2','MRPL2','UQCRQ','MRPL15','TMEM70','MRPS16','MRPL43','MRPS26','NDUFB3','BCS1L','UQCC3','TTC19','POLRMT'),
    'aur' = c('MRPS21','NDUFAF3','ATP5MF','UQCRQ','MRPL43','MRPL51','BCS1L','UQCC3','MRPL10','CCDC56','CARS2'),
    'fas_2' = c('NDUFAF3','MRPL18','RARS2','MRPL2','NDUFAF2','MRPS36','MRPS27','NDUFA2','CD14','NDUFAF6','CYC1','BCS1L','MRPL53','COX7C','MRPL49','NARS2','TMEM126B','MRPL10','CCDC56','NDUFA13','NDUFB10'),
    'thi_2' = c('MRPL18','MRPL2','NDUFA2','CD14','MRPL49','LOC722212','MRPL10','CCDC56','NDUFB10'),
    'ass_2' = c('ATP5MF','MRPL36','NDUFS6','UQCC3','MRPL12','CARS2','POLRMT','NDUFA13','MRPS34','NDUFB10')
)

