library(ggrepel)
library(ggplot2)
library(jtools) # for summ()
library(interactions)
library(sjPlot) 
library(ggeffects) 
library(dbplyr)
library(pointr)
library(emmeans)
library(gridExtra)
library(patchwork)
library(gtable)
library(grid)
library(ggpubr)

drawing_predict <- function(lmodel,tag_here, x_limit){
    to_draw <- paste0("num_genes [1:", x_limit,"]")
    temp <- ggpredict(lmodel,  c(to_draw,"containsNinteractgenez"), back.transform = F)
    y_max = ceiling(max(lmodel$model$`log(length)`[lmodel$model$num_genes == x_limit]))
    y_min = floor(min(lmodel$model$`log(length)`[lmodel$model$num_genes == 1]))
    attr(temp, "rawdata")$response <- lmodel$model$`log(length)`
    attr(temp, "rawdata") <- attr(temp, "rawdata")[order(attr(temp, "rawdata")$group, decreasing = F), ]
    
    plot(temp, show.y.title = F, show.title = FALSE,
         show.x.title = F,
         show.legend = F,
         add.data = TRUE,
         colors = c("blue","red"),
         limit.range=F,
    )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   scale_x_continuous(name="", limits=c(0,x_limit), breaks = 0:x_limit) +
   scale_y_continuous(name="", limits=c(y_min,y_max)) +
    theme(axis.text.x = element_text(size=14))+
    theme(axis.title.x = element_text(size=16))+
    theme(axis.text.y = element_text(size=14))+
    theme(axis.title.y = element_text(size=16)) + 
    theme(plot.tag = element_text(size = 16)) +
    labs(tag = tag_here)
}

drawing_predict <- function(lmodel,tag_here, x_limit){
    ggpredict(lmodel,  c("num_genes","containsNinteractgenez"), back.transform = F) %>%
        plot(show.y.title = FALSE, show.title = FALSE,
             show.x.title = FALSE,
             show.legend = F,
             colors = c("blue","red"),
             add.data = TRUE, #### this adds points to the plot  ####
             limit.range=F
        )+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        scale_x_continuous(name="", limits=c(0,x_limit)) +
#        scale_y_continuous(name=" ", limits=c(0,3), breaks = c(0, 50, 100, 150))+
        theme(axis.text.x = element_text(size=14))+
        theme(axis.title.x = element_text(size=16))+
        theme(axis.text.y = element_text(size=14))+
        theme(axis.title.y = element_text(size=16)) +
        labs(tag = tag_here)
}    


justify <- function(x, hjust="center", vjust="top", draw=FALSE){
    w <- sum(x$widths)
    h <- sum(x$heights)
    xj <- switch(hjust,
                 center = 0.5,
                 left = 0.5*w,
                 right=unit(1,"npc") - 0.5*w)
    yj <- switch(vjust,
                 center = 0.5,
                 bottom = 0.5*h,
                 top=unit(1,"npc") - 0.8*h)
    x$vp <- viewport(x=xj, y=yj)
    if(draw) grid.draw(x)
    return(x)
}


setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/roh_results/")

theme_set(theme_ggeffects())

species_names <- c('M.lasiotis','M.littoralis','M.brevicaudus','M.mulatta','M.tcheliensis')
species_names <- c('Blue','Red','Purple')
species_names <- c('JC_UoC','RWNK_WNPRC','LPVGTI_OHCU','EV_NEPRC','ZJ_YNPRC','BF_ONPRC','MK_TNPRC','SKDS_CNPRC','JH_CPRC','RWDO_WNPRC')
species_names <- c("arctoides","assamensis","thibetana","mulata","fascicularis")
species_names <- c("aureus","fascicularis","thibetana","assamensis")
species_names <- c('orange','brown','red')
root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/chinese_mito"

for (species in species_names){
    # read.table(file.path(root,'ROH_Blue_roh_density.out'), header=T)
    assign(species,read.table(file.path(root,paste0("roh_",species,".density.out")), header = T)) 
    ptr("now_data", species)
    now_data <- now_data[now_data$containsgenes == 1,]
    now_data$length <- now_data$length/100000
    now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 
    
    # species_filter <- paste0(species,"_allgenez")
    # assign(species_filter,now_data[(now_data$containsgenes == 1),])
    # ptr("now_data_filtered", species_filter)
    # 
    # thresh_hold <- quantile(now_data_filtered$num_genes, probs = seq(0.99,0.99))
    # now_data_filtered <- now_data_filtered[now_data_filtered$num_genes < thresh_hold,]
    
    assign(paste0(species,"_mod"),lm(log(length) ~ containsNinteractgenez * num_genes, data= now_data))
}

vec <- paste0("(",letters,")")
cnt = 0
for (species in species_names) {
    cnt = cnt +1
    a = eval(parse(text=paste(species,"mod",sep="_")))
    x_limit = quantile(a$model$num_genes, probs = seq(0.9, 0.9))
    assign(paste(species,'plot',sep="_"), drawing_predict(eval(parse(text=paste(species,"mod",sep="_"))),vec[cnt], x_limit))
}

writefile = "quantile_prop.txt"
sink(writefile)
for (species in species_names) {
    print(paste0("# ", species))
    print(quantile(eval(parse(text=paste(species)))$num_genes[eval(parse(text=paste(species)))$num_genes > 0], probs = seq(0, 1, 0.1)))
}
sink()

print(thibetana_plot)
summary(thibetana_mod)

bottom_text = text_grob("Number of genes",size = 20,just = 'center')
bottom_text$vp <- viewport(x=0.5, y=1)

left_text = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 1,size = 20)
left_text$vp <- viewport(x=0.2, y=0.2)

# 10, 6
png("Predicted_values_ROH_aureus.png", width = 168, height = 200, units='mm', res = 300)    
grid.arrange(
    arrangeGrob(
        grobs = list(aureus_plot,fascicularis_plot,thibetana_plot,assamensis_plot),
        ncol = 2
    ),
    left = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 0.5,hjust = 0.5,size = 20),
    bottom = bottom_text
)
dev.off()

# 10, 5
png("Predicted_values_ROH_arctoides.png", width = 250, height = 200, units='mm', res = 300)    
grid.arrange(
    arrangeGrob(
        grobs = list(arctoides_plot,assamensis_plot,thibetana_plot,mulata_plot,fascicularis_plot),
        ncol = 3
    ),
    left = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 0.5,hjust = 0.5,size = 20),
    bottom = bottom_text
)
dev.off()

# 10, 3
png("Predicted_values_ROH_chinese.png", width = 250, height = 100, units='mm', res = 300)    
grid.arrange(
    arrangeGrob(
        grobs = list(Blue_plot,Red_plot,Purple_plot),
        ncol = 3
    ),
    left = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 0.5,hjust = 0.5,size = 20),
    bottom = bottom_text
)
dev.off()

png("Predicted_values_ROH_indian.png", width = 250, height = 100, units='mm', res = 300)    
grid.arrange(
    arrangeGrob(
        grobs = list(orange_plot,brown_plot, red_plot),
        ncol = 3
    ),
    left = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 0.5,hjust = 0.5,size = 20),
    bottom = bottom_text
)
dev.off()
