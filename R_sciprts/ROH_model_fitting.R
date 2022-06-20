library(ggrepel)
library(ggplot2)
library(jtools) # for summ()
library(interactions)
library(sjPlot) 
library(ggeffects) 
library(dbplyr)
library(pointr)
library(emmeans)

drawing_interation <- function(lmodel){
  ggpredict(lmodel,  c("num_genes","containsNinteractgenez")) %>%
    plot(show.y.title = F, show.title = FALSE,
         show.x.title = F,
         show.legend = F,
         colors = c("blue","red"),
         limit.range=F
    )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    # scale_x_continuous(name="", limits=c(0,50)) +
    
    theme(axis.text.x = element_text(size=14))+
    theme(axis.title.x = element_text(size=16))+
    theme(axis.text.y = element_text(size=14))+
    theme(axis.title.y = element_text(size=16))
}

drawing_distribution <- function(dataf){
  ggplot(dataf) + 
    geom_density(aes(x=num_genes, colour=group_f))+
    stat_density(aes(x=num_genes, colour=group_f),
                 geom="line",position="identity")+
    scale_color_manual(values = c("blue", "red"))+
    theme_bw()
}

drawing_relation <- function(dataf,dataf1,name){
  plot(x = dataf[dataf$containsNinteractgenez == 0, ]$num_genes, 
       y = dataf[dataf$containsNinteractgenez == 0, ]$length, 
       col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25), pch = 19,
       xlab = "Gene numbers", ylab = "ROH length", main= species)
  abline(a = coef(dataf1)[1], b = coef(dataf1)[3], col = "blue", pch = 19, lwd = 2) +
    points(x = dataf[dataf$containsNinteractgenez == 1, ]$num_genes, 
           y = dataf[dataf$containsNinteractgenez == 1, ]$length, 
           col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19) +
    abline(a = coef(dataf1)[1] + coef(dataf1)[2], b = coef(dataf1)[3] + coef(dataf1)[4], 
           col = "red", lwd = 2)
}

setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/roh_results/indian_mito")

theme_set(theme_ggeffects())

species_names <- c('M.lasiotis','M.littoralis','M.brevicaudus','M.mulatta','M.tcheliensis')
species_names <- c('Blue','Red','Purple')
species_names <- c('JC_UoC','RWNK_WNPRC','LPVGTI_OHCU','EV_NEPRC','ZJ_YNPRC','BF_ONPRC','MK_TNPRC','SKDS_CNPRC','JH_CPRC','RWDO_WNPRC')
species_names <- c("arctoides","assamensis","thibetana","mulata","fascicularis")
species_names <- c("aureus","fascicularis","thibetana","assamensis")
species_names <- c('orange','brown','red')
root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/indian_mito"

# root <- "G:/Now_work/2021_05_18.ben_/09.re_analysis/roh/aureus"
# species_names <- c("Chinese","Indian","Red")

write_file <- file("lm_roh_new.md", encoding = "UTF-8")
sink(write_file,append = T)


# paste("roh_",species,"density.out", sep=".")

# Write output of linear models
for (species in species_names){
  # read.table(file.path(root,'ROH_Blue_roh_density.out'), header=T)
  assign(species,read.table(file.path(root,paste0("roh_",species,".density.out")), header = T)) 
  ptr("now_data", species)
  now_data$length <- now_data$length/100000
  now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 

  species_filter <- paste0(species,"_allgenez")
  assign(species_filter,now_data[(now_data$containsgenes == 1),])
  ptr("now_data_filtered", species_filter)

  assign(paste0(species,"_genes_or_not"),lm(log(length) ~ containsgenes, data=now_data))
  assign(paste0(species,"_number_of_genes_all_ROHs"),lm(log(length) ~ num_genes, data=now_data))
  assign(paste0(species,"_number_of_genes"),lm(log(length) ~ num_genes, data=now_data_filtered))
  assign(paste0(species,"_Ninteract"),lm(log(length) ~ containsNinteractgenez, data=now_data_filtered))
  assign(paste0(species,"_mod"),lm(log(length) ~ containsNinteractgenez*num_genes, data=now_data_filtered))
  
  cat(paste("##",species))
  cat("\n")
  cat("### overall relation")
  cat("\n")
  cat('#### containsgenes')
  cat("\n")
  cat('##### results')
  cat("\n")
  print(summary(eval(parse(text=paste(species,"genes_or_not",sep="_")))))
  cat("\n")
  cat('#### density')
  cat("\n")
  cat('##### results')
  cat("\n")
  print(summary(eval(parse(text=paste(species,"number_of_genes_all_ROHs",sep="_")))))
  cat("\n")
  cat('### only genic ROHs')
  cat("\n")
  cat('#### density')
  cat("\n")
  cat('##### results')
  cat("\n")
  print(summary(eval(parse(text=paste(species,"number_of_genes",sep="_")))))
  cat("\n")
  cat('#### N-interact')
  cat("\n")
  cat('##### results')
  cat("\n")
  print(summary(eval(parse(text=paste(species,"Ninteract",sep="_")))))
  cat("\n")
  cat('#### N-interact * density')
  cat("\n")
  cat('##### results')
  cat("\n")
  print(summary(eval(parse(text=paste(species,"mod",sep="_")))))
  cat("\n")
  cat("\n")
  
  # to_do_list = c(
  #   paste0(species,"_genes_or_not"),paste0(species,"_number_of_genes"),paste0(species,"_number_of_genes_all_ROHs"),
  #   paste0(species,"_Ninteract"),paste0(species,"_mod")
  # )
  # 
  # save(list=to_do_list,file = paste0(species,".rs"))
  # rm(list=to_do_list)
  # 
}
sink()

# Release memory
for (species in species_names) {
  load(paste0(species,".rs"))
  to_do_list = c(
    paste0(species,"_genes_or_not"),paste0(species,"_number_of_genes"),paste0(species,"_number_of_genes_all_ROHs"),
    paste0(species,"_Ninteract")
  )
  rm(list=to_do_list)
}

# species <- species_names[4]
# load(paste0(species,".rs"))



###

for (species in species_names) {
  assign(paste(species,'plot',sep="_"), drawing_interation(eval(parse(text=paste(species,"mod",sep="_")))))
}

vec <- paste0("(",letters,")")

png("Predicted_values_ROH_ln.png", width = 250, height = 200, units='mm', res = 300)    
  # plot_grid(list(M.lasiotis_plot, M.littoralis_plot, M.brevicaudus_plot,
  #                M.mulatta_plot, M.tcheliensis_plot
  # ),
 plot_grid(list(orange_plot,brown_plot, red_plot),
  # plot_grid(list(JC_UoC_plot,RWNK_WNPRC_plot,LPVGTI_OHCU_plot,EV_NEPRC_plot,ZJ_YNPRC_plot,BF_ONPRC_plot,MK_TNPRC_plot,SKDS_CNPRC_plot,JH_CPRC_plot,RWDO_WNPRC_plot),
#  plot_grid(list(arctoides_plot,assamensis_plot,thibetana_plot,mulata_plot,fascicularis_plot),
#    sjPlot::plot_grid(list(aureus_plot,fascicularis_plot,thibetana_plot,assamensis_plot),
#  plot_grid(list(Blue_plot,Red_plot,Purple_plot),
  tags = vec[1:length(species_names)],
  margin = c(0.1, 0.1, 0.1, 0.1))+
    theme_bw() + ggplot2::xlab("gene")
dev.off()

###

#all_data_allgenez <- rbind(dataf,M.lasiotis_allgenez,M.littoralis_allgenez,
#                           M.mulatta_allgenez,M.tcheliensis_allgenez)
all_data_allgenez <- rbind(M.lasiotis_allgenez,M.littoralis_allgenez,M.mulatta_allgenez,M.tcheliensis_allgenez,M.brevicaudus_allgenez)
all_data_allgenez <- rbind(JC_UoC_allgenez,RWNK_WNPRC_allgenez,LPVGTI_OHCU_allgenez,EV_NEPRC_allgenez,ZJ_YNPRC_allgenez,BF_ONPRC_allgenez,MK_TNPRC_allgenez,SKDS_CNPRC_allgenez,JH_CPRC_allgenez,RWDO_WNPRC_allgenez)
all_data_allgenez <- rbind(arctoides_allgenez,assamensis_allgenez,thibetana_allgenez,mulata_allgenez,fascicularis_allgenez)
all_data_allgenez <- rbind(aureus_allgenez,fascicularis_allgenez,thibetana_allgenez,assamensis_allgenez)
all_data_allgenez <- rbind(Blue_allgenez,Purple_allgenez,Red_allgenez)
all_data_allgenez <- rbind(green_allgenez,brown_allgenez,special_allgenez)

png("histo_of_all_data_density.png", width = 250, height = 200, units='mm', res = 300) 
hist(all_data_allgenez$gene_density_on_ROHs,xlim=c(0,6),warn.unused=F) #breaks=seq(-1,2000,by=1)
dev.off()
breaks=seq(-1,152,by=1)
# calculate the proportion of ROHs that have a gene density above some value
length(all_data_allgenez$gene_density_on_ROHs[all_data_allgenez$gene_density_on_ROHs > 20])/length(all_data_allgenez$gene_density_on_ROHs)


# calculate the proportion of ROHs that have a gene density above some value
length(all_data_allgenez$gene_density_on_ROHs[all_data_allgenez$gene_density_on_ROHs < 5])/length(all_data_allgenez$gene_density_on_ROHs)

all_data <- data.frame()
for (species in species_names) {
  print(species)
  print(eval(parse(text=paste0("max(",species,"$num_ninteractgenes)"))))
  overall_max <- eval(parse(text=paste0("max(",species,"$num_ninteractgenes)")))
  print('contains')
  print(eval(parse(text=paste0("mean(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==1], na.rm=T)"))))
  contain_mean <- eval(parse(text=paste0("mean(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==1], na.rm=T)")))
  print(eval(parse(text=paste0("max(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==1], na.rm=T)"))))
  contain_max <- eval(parse(text=paste0("max(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==1], na.rm=T)")))
  print('not_contains')
  print(eval(parse(text=paste0("mean(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==0], na.rm=T)"))))
  nc_mean <- eval(parse(text=paste0("mean(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==0], na.rm=T)")))
  print(eval(parse(text=paste0("max(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==0], na.rm=T)"))))
  nc_max <- eval(parse(text=paste0("max(",species,"_allgenez$gene_density_on_ROHs[",species,"_allgenez$containsNinteractgenez ==0], na.rm=T)")))
  temp_frame <-data.frame(overall_max,contain_mean,contain_max,nc_mean,nc_max,row.names = species)
  all_data <- rbind(all_data, temp_frame)
}

write.csv(all_data, quote=F,file = 'data_summary.csv')


for (species in species_names) {
  cmd <- paste0(species,"_allgenez$group_f = factor(",species,"_allgenez$containsNinteractgenez, levels=c('0','1'), ordered = T)")
  eval(parse(text=cmd))
  assign(paste0(species,"_distribution_plot"),drawing_distribution(eval(parse(text=paste0(species,"_allgenez")))))
}

png("Distributions.png", width = 250, height = 200, units='mm', res = 1000)    
#plot_grid(list(M.lasiotis_distribution_plot, M.littoralis_distribution_plot, M.brevicaudus_distribution_plot,
#               M.mulatta_distribution_plot, M.tcheliensis_distribution_plot),
#plot_grid(list(arctoides_distribution_plot,assamensis_distribution_plot,thibetana_distribution_plot,mulata_distribution_plot,fascicularis_distribution_plot),
plot_grid(list(aureus_distribution_plot,fascicularis_distribution_plot,thibetana_distribution_plot,assamensis_distribution_plot),
#  plot_grid(list(Blue_distribution_plot,Red_distribution_plot,Purple_distribution_plot),
#  plot_grid(list(green_distribution_plot,brown_distribution_plot, special_distribution_plot),
tags = species_names,
margin = c(0.1, 0.1, 0.1, 0.1))+
  theme_bw()
dev.off()

png("Distributions.png", width = 500, height = 200, units='mm', res = 1000)  
plot_grid(list(JC_UoC_distribution_plot,RWNK_WNPRC_distribution_plot,LPVGTI_OHCU_distribution_plot,EV_NEPRC_distribution_plot,ZJ_YNPRC_distribution_plot,BF_ONPRC_distribution_plot,MK_TNPRC_distribution_plot,SKDS_CNPRC_distribution_plot,JH_CPRC_distribution_plot,RWDO_WNPRC_distribution_plot
),
  tags = species_names,
  margin = c(0.1, 0.1, 0.1, 0.1)
)+
  theme_bw()
dev.off()

plot_grid(list(Blue_distribution_plot, Red_distribution_plot, Purple_distribution_plot),
          tags = species_names,
          margin = c(0.1, 0.1, 0.1, 0.1)          
  ) +
  theme_bw()
dev.off()

for (species in species_names) {
  emtrends(M.mulatta_mod, pairwise ~ containsNinteractgenez, var="gene_density_on_ROHs")
  mylist <- list(gene_density_on_ROHs=seq(0,80,by=10))
  emmip(M.mulatta_mod, containsNinteractgenez~gene_density_on_ROHs, at=mylist, CIs=TRUE)
}



png("Density_length.png", width = 500, height = 200, units='mm', res = 1000)  
par(mfrow=c(2,2))
for (species in species_names) {
  drawing_relation(eval(parse(text=paste0(species,"_allgenez"))),eval(parse(text=paste0(species,"_mod"))),species)
}
dev.off()

plot(c(1,2,3))
par(xpd=TRUE)
plot_colour <- c("red","blue")
text_here <- c("with N-interact","wihtout N-interact")
legend(x="bottom",legend = text_here,col=plot_colour,lwd=5, cex=1)
par(xpd=F)

par(mfrow=c(1,1))
################################################

# my_data_only_genez <- mau[mau$containsgenes == 1,] 

# # What are the names of the Ninteract genes with no polymorphism (FW_H = 'NA')
# my_data_only_genez$Ninteract_acronym[((is.na(my_data_only_genez$FW_H)) &
#                                         (my_data_only_genez$containsNinteractgenez == 1))]
# # ACAD9

# # explore relationship between FW_H and number of genes in Ninteract windows
# my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
# my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
# my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]

# my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
# #dim(my_data_only_genez)
# #head(my_data_only_genez)

# # calculate a lm for all data (because there was not a significant interaction term)
# mod <- lm(FW_H ~ number_of_genes, data=my_data_only_genez)
# # get the fitted values (y = mx+b)
# fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
# #cbind fitted to data
# my_data_only_genez <- cbind(my_data_only_genez,fitted)

# # calculate cooks d for all data
# cooksd <- cooks.distance(mod)
# #cbind cooksd to data
# my_data_only_genez <- cbind(my_data_only_genez,cooksd)
# # make a column to specify whether a gene is an Ninteract gene or not
# my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
# my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
# # make a column that specifies whether cooksd suggests an outlier
# # but only for Ninteract genes
# my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
#                            (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
# # now color the ones that are below the fitted line blue
# my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
#                            (my_data_only_genez$FW_H < my_data_only_genez$fitted) &
#                            (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"

# # get some numbers
# # What is the expected proportion of upper outliers
# nrow(my_data_only_genez[(cooksd >=4*mean(cooksd, na.rm=T)) & # cooksD is high
#                           (my_data_only_genez$containsNinteractgenez == 0) & # the gene is not an Ninteract gene
#                           (my_data_only_genez$FW_H < my_data_only_genez$fitted),])/ # the value is an upper outlier
#   nrow(my_data_only_genez[(my_data_only_genez$containsNinteractgenez == 0),])
# #  0.02516556
# # What is the observed proportion of upper outliers
# nrow(my_data_only_genez[(cooksd >=4*mean(cooksd, na.rm=T)) & # cooksD is high
#                           (my_data_only_genez$containsNinteractgenez == 1) & # the gene is  an Ninteract gene
#                           (my_data_only_genez$FW_H < my_data_only_genez$fitted),])/ # the value is an upper outlier
#   nrow(my_data_only_genez[(my_data_only_genez$containsNinteractgenez == 1),])
# #  0.07881773

# # what are the names of the Ninteract outlier windows?
# my_data_only_genez$Ninteract_acronym[(cooksd >=4*mean(cooksd, na.rm=T)) &
#                                        (my_data_only_genez$containsNinteractgenez == 1) &
#                                        (my_data_only_genez$FW_H < my_data_only_genez$fitted)]

# #  MRPL55        MRPL53        MRPL30        NDUFAF3       UQCRQ         MRPL18        ATP5J2        C11orf83     
# # MRPL40        MRPS34        EARS2,NDUFAB1 COA3          MRPL38        NDUFS7        COX6B1        ATP5SL   


# # make the color column into an ordered factor
# #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
# #                                   ordered = is.ordered(my_data_only_genez$color))
# # on now plot the data with the color representing outliers for N_interact genes
# mau_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = FW_H, color=color, fill=color)) +
#   geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
#   #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
#   #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
#   geom_point(data = subset(my_data_only_genez, color == "gray"),
#              aes(x = number_of_genes, y = FW_H, alpha = alpha), color = 'gray') +
#   geom_point(data = subset(my_data_only_genez, color == 'pink'),
#              aes(x = number_of_genes, y = FW_H, alpha = alpha), color = 'pink') +
#   geom_point(data = subset(my_data_only_genez, color == 'red'),
#              aes(x = number_of_genes, y = FW_H, alpha = alpha),color = 'red') +
#   geom_point(data = subset(my_data_only_genez, color == 'blue'),
#              aes(x = number_of_genes, y = FW_H, alpha = alpha),color = 'red') +
#   #      geom_text_repel( data = my_data_only_genez,
#   #                       #mapping = aes(label = Ninteract_acronym),
#   #                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
#   #                       force_pull = 0,
#   #                       force = 13,
#   #                       nudge_y = 0.1, nudge_x = 15,
#   #                       color = "black",
#   #                       size = 2.5,
#   #                       box.padding = 0.5, 
#   #                       #point.padding = 0.5,
#   #                       direction     = "y",
# #                       max.overlaps = Inf,
# #                       hjust = 0, #angle = 45, #segment.curvature = -0.05,
# #                       segment.size = 0.25,
# #                       segment.color = 'grey50'
# #      ) +
# xlim(0,16) + scale_y_continuous(limits = c(-8,6), breaks = c(-6.0,-3.0,0,3.0,6.0)) +
#   labs(x = element_blank(), y="H", tag = "mau") +
#   theme_classic(base_size=16) + theme(legend.position = "none") +
#   theme(axis.text.x=element_blank()) + #theme(axis.text.y=element_blank()) +
#   theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

# summary(glm(length ~ containsNinteractgenez, data=M.brevicaudus))
