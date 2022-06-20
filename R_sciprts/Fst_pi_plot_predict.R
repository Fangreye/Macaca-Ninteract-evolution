library(ggrepel)
library(ggplot2)
library(jtools) # for summ()
library(interactions)
library(sjPlot) 
library(ggeffects) 
library(dbplyr)
library(pointr)
library(emmeans)
library(egg)
library(gridExtra)
library(patchwork)
library(gtable)
library(grid)

root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/'
setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/roh_results")
vec <- paste0("(",letters,")")

######### For Fst

root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/new_fst_out/'
setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/new_fst_out")
vec <- paste0("(",letters,")")    

#species <- c("aureus_new","chinese_mito","chinese_org","fas","indian_mito","indian_org")
species <- c("aureus_new","chinese_mito","fas","indian_mito")
pairses_list <- list(
    #c('aureus_fascicularis','thibetana_assamensis_aureus','thibetana_assamensis_fascicularis'),
    c('aureus_assamensis','aureus_thibetana','aureus_fascicularis','fascicularis_assamensis','fascicularis_thibetana','assamensis_thibetana'),
    c('Blue_Purple','Blue_Red','Purple_Red'),
#    c('M.brevicaudus_M.lasiotis','M.brevicaudus_M.mulatta','M.lasiotis_M.mulatta',
      # 'M.littoralis_M.brevicaudus','M.littoralis_M.lasiotis','M.littoralis_M.mulatta',
      # 'M.tcheliensis_M.brevicaudus','M.tcheliensis_M.lasiotis','M.tcheliensis_M.littoralis',
      # 'M.tcheliensis_M.mulatta'),
    c ("arctoides_assamensis","arctoides_thibetana","arctoides_mulata","arctoides_fascicularis","assamensis_thibetana","assamensis_mulata","assamensis_fascicularis","thibetana_mulata","thibetana_fascicularis","mulata_fascicularis"),
    #  c('brown_green','brown_special','orange_special'),
    c('brown_orange','brown_special','orange_special')
#    c('BF_ONPRC_JH_CPRC','BF_ONPRC_MK_TNPRC','BF_ONPRC_RWDO_WNPRC','BF_ONPRC_SKDS_CNPRC','EV_NEPRC_BF_ONPRC','EV_NEPRC_JH_CPRC','EV_NEPRC_MK_TNPRC','EV_NEPRC_RWDO_WNPRC','EV_NEPRC_SKDS_CNPRC','EV_NEPRC_ZJ_YNPRC','JC_UoC_BF_ONPRC','JC_UoC_EV_NEPRC','JC_UoC_JH_CPRC','JC_UoC_LPVGTI_OHCU','JC_UoC_MK_TNPRC','JC_UoC_RWDO_WNPRC','JC_UoC_RWNK_WNPRC','JC_UoC_SKDS_CNPRC','JC_UoC_ZJ_YNPRC','JH_CPRC_RWDO_WNPRC','LPVGTI_OHCU_BF_ONPRC','LPVGTI_OHCU_EV_NEPRC','LPVGTI_OHCU_JH_CPRC','LPVGTI_OHCU_MK_TNPRC','LPVGTI_OHCU_RWDO_WNPRC','LPVGTI_OHCU_SKDS_CNPRC','LPVGTI_OHCU_ZJ_YNPRC','MK_TNPRC_JH_CPRC','MK_TNPRC_RWDO_WNPRC','MK_TNPRC_SKDS_CNPRC','RWNK_WNPRC_BF_ONPRC','RWNK_WNPRC_EV_NEPRC','RWNK_WNPRC_JH_CPRC','RWNK_WNPRC_LPVGTI_OHCU','RWNK_WNPRC_MK_TNPRC','RWNK_WNPRC_RWDO_WNPRC','RWNK_WNPRC_SKDS_CNPRC','RWNK_WNPRC_ZJ_YNPRC','SKDS_CNPRC_JH_CPRC','SKDS_CNPRC_RWDO_WNPRC','ZJ_YNPRC_BF_ONPRC','ZJ_YNPRC_JH_CPRC','ZJ_YNPRC_MK_TNPRC','ZJ_YNPRC_RWDO_WNPRC','ZJ_YNPRC_SKDS_CNPRC')
    
)

species <- c("aureus_new")
pairses_list <- list(
    #c('aureus_fascicularis','thibetana_assamensis_aureus','thibetana_assamensis_fascicularis'),
    c('assamensis_DRR219369','assamensis_DRR219370','assamensis_fascicularis',
      'DRR219369_DRR219370','DRR219369_fascicularis','DRR219370_fascicularis',
      'thibetana_assamensis','thibetana_DRR219369','thibetana_DRR219370',
      'thibetana_fascicularis')
)

for (size in c("30k")) {
    for (i in seq(1,4)) {
        speciy <- species[i]
        pairses <- pairses_list[[i]]
        dat_list <- rep("",length(pairses))
        
        long_list = "rbind("
        for (pair in pairses) {
            
            point_h <- paste(pair,'roh',sep='.')
            
            filpath <- file.path( root,'fst_out',size,speciy,paste0("fst_",pair,'.density.out'))
            assign(point_h,read.table(filpath, header = T))
            
            cmd <- paste0(point_h," <- data.frame(",point_h," , \"species\" = rep(\"",pair,"\", dim(",point_h,")[1]))")
            eval(parse(text = cmd))
            
            long_list = paste0(long_list,point_h,",")
            
        }
        
        long_list <- paste0(substr(long_list, 1,nchar(long_list)-1),")")
        
        sel <- eval(parse(text = long_list))
        
        sel$group <- "NA"
        sel$group[which((sel$containsgenes == 0)&(sel$containsNinteractgenez == 0))] <- "No genes"
        sel$group[which((sel$containsgenes == 1)&(sel$containsNinteractgenez == 0))] <- "Other genes"
        sel$group[which((sel$containsgenes == 1)&(sel$containsNinteractgenez == 1))] <- "Ninteract genes"
        
        sel$group_f = factor(sel$group, levels=c('No genes','Other genes','Ninteract genes'), ordered = T)
        sel$species_f = factor(sel$species, levels=pairses, labels = vec[1:length(pairses)], ordered = T)
        
        paa <- ggplot(sel) + 
            geom_density(aes(x=Fst, colour=group_f),show_guide=FALSE)+
            stat_density(aes(x=Fst, colour=group_f),
                         geom="line",position="identity")+
            scale_color_manual(values = c("black", "blue", "red")) +
            facet_wrap(~species_f, scales = "free", ncol=3) +
            theme_bw() +
            # remove legend key border color & background
            theme(legend.key=element_blank()) +
            theme(legend.box.background = element_blank())+
            theme(legend.title = element_blank()) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            theme(strip.background = element_blank()) +
            labs(x = expression(italic("F"["ST"])), y = "Density") + 
            theme(
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 22),
                axis.title.x = element_text(size = 28),
                axis.title.y = element_text(size = 28),
                strip.text = element_text(size = 26),
                legend.text = element_text(size = 26)) +
            theme(legend.position = c(0.09, 1 - 0.2 /ceiling(length(pairses)/3))) +
            theme(legend.background=element_rect(fill = alpha("white", 0))) + 
            theme(panel.spacing.x = unit(8, "mm"))
#        print(paa)
#        dev.off()
        
    }    
    
}

######### For pi

#species <- c("aureus_new","chinese_mito","chinese_org","fas","indian_mito","indian_org")
species <- c("aureus_new","chinese_mito","fas","indian_mito")
pairses_list <- list(
    #c('aureus','fascicularis','thibetana_assamensis'),
    c('aureus','fascicularis','thibetana','assamensis'),
    c('Blue','Red','Purple'),
#    c('M.brevicaudus','M.lasiotis','M.littoralis','M.mulatta','M.tcheliensis'),
    c("arctoides","assamensis","thibetana","mulata","fascicularis"),
    c('orange','brown','red')
#    c('JC_UoC','RWNK_WNPRC','LPVGTI_OHCU','EV_NEPRC','ZJ_YNPRC','BF_ONPRC','MK_TNPRC','SKDS_CNPRC','JH_CPRC','RWDO_WNPRC')
)

for (size in c("30k")) {
    for (i in seq(1,4)) {
        speciy <- species[i]
        pairses <- pairses_list[[i]]
        dat_list <- rep("",length(pairses))
        
        long_list = "rbind("
        for (pair in pairses) {
            
            point_h <- paste(pair,'roh',sep='.')
            
            filpath <- file.path( root,'pi_out',size,speciy,paste0("pi_",pair,'.density.out'))
            #filpath <- file.path( root,'fst_out',size,speciy,paste0("fst_",pair,'.density.out'))
            assign(point_h,read.table(filpath, header = T))
            
            cmd <- paste0(point_h," <- data.frame(",point_h," , \"species\" = rep(\"",pair,"\", dim(",point_h,")[1]))")
            eval(parse(text = cmd))
            
            long_list = paste0(long_list,point_h,",")
            
        }
        
        long_list <- paste0(substr(long_list, 1,nchar(long_list)-1),")")
        
        sel <- eval(parse(text = long_list))
        
        sel$group <- "NA"
        sel$group[which((sel$containsgenes == 0)&(sel$containsNinteractgenez == 0))] <- "No genes"
        sel$group[which((sel$containsgenes == 1)&(sel$containsNinteractgenez == 0))] <- "Other genes"
        sel$group[which((sel$containsgenes == 1)&(sel$containsNinteractgenez == 1))] <- "Ninteract genes"
        
        sel$group_f = factor(sel$group, levels=c('No genes','Other genes','Ninteract genes'), ordered = T)
        sel$species_f = factor(sel$species, levels=pairses, labels = vec[1:length(pairses)], ordered = T)
        
        png(paste("Pi_density_allNinteract_",speciy,"_",size,".png",sep=""),
            width = 500, height = ceiling(118*length(pairses)/3), units='mm', res = 100)
        paa <- ggplot(sel) + 
            geom_density(aes(x=pi, colour=group_f),show_guide=FALSE)+
            stat_density(aes(x=pi, colour=group_f),
                         geom="line",position="identity")+
            scale_color_manual(values = c("black", "blue", "red")) +
            facet_wrap(~species_f, scales = "free", ncol=3) +
            theme_bw() +
            # remove legend key border color & background
            theme(legend.key=element_blank()) +
            theme(legend.box.background = element_blank())+
            theme(legend.title = element_blank()) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            theme(strip.background = element_blank()) +
            labs(x = expression(italic("\u03c0")), y = "Density") + 
            theme(
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 22),
                axis.title.x = element_text(size = 28,family="serif"),
                axis.title.y = element_text(size = 28),
                strip.text = element_text(size = 26),
                legend.text = element_text(size = 26)) +
            theme(legend.position = c(0.09, 1 - 0.2 /ceiling(length(pairses)/3))) +
            theme(legend.background=element_rect(fill = alpha("white", 0))) + 
            theme(panel.spacing.x = unit(8, "mm"))
        print(paa)
        dev.off()
    }    
    
}
