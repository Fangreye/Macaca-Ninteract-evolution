library(genomic.autocorr)

setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/fst_out")
root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/'

word <- c("100k","30k")
species <- c("aureus_new","chinese_mito","chinese_org","fas","indian_mito","indian_org")

# For Fst
pairses_list <- list(
    #c('aureus_fascicularis','thibetana_assamensis_aureus','thibetana_assamensis_fascicularis'),
    c('aureus_assamensis','aureus_thibetana','aureus_fascicularis','fascicularis_assamensis','fascicularis_thibetana','assamensis_thibetana'),
    c('Blue_Purple','Blue_Red','Purple_Red'),
    c('M.brevicaudus_M.lasiotis','M.brevicaudus_M.mulatta','M.lasiotis_M.mulatta',
      'M.littoralis_M.brevicaudus','M.littoralis_M.lasiotis','M.littoralis_M.mulatta',
      'M.tcheliensis_M.brevicaudus','M.tcheliensis_M.lasiotis','M.tcheliensis_M.littoralis',
      'M.tcheliensis_M.mulatta'),
    c ("arctoides_assamensis","arctoides_thibetana","arctoides_mulata","arctoides_fascicularis","assamensis_thibetana","assamensis_mulata","assamensis_fascicularis","thibetana_mulata","thibetana_fascicularis","mulata_fascicularis"),
    #  c('brown_green','brown_special','orange_special'),
    c('brown_orange','brown_special','orange_special'),
    c('BF_ONPRC_JH_CPRC','BF_ONPRC_MK_TNPRC','BF_ONPRC_RWDO_WNPRC','BF_ONPRC_SKDS_CNPRC','EV_NEPRC_BF_ONPRC','EV_NEPRC_JH_CPRC','EV_NEPRC_MK_TNPRC','EV_NEPRC_RWDO_WNPRC','EV_NEPRC_SKDS_CNPRC','EV_NEPRC_ZJ_YNPRC','JC_UoC_BF_ONPRC','JC_UoC_EV_NEPRC','JC_UoC_JH_CPRC','JC_UoC_LPVGTI_OHCU','JC_UoC_MK_TNPRC','JC_UoC_RWDO_WNPRC','JC_UoC_RWNK_WNPRC','JC_UoC_SKDS_CNPRC','JC_UoC_ZJ_YNPRC','JH_CPRC_RWDO_WNPRC','LPVGTI_OHCU_BF_ONPRC','LPVGTI_OHCU_EV_NEPRC','LPVGTI_OHCU_JH_CPRC','LPVGTI_OHCU_MK_TNPRC','LPVGTI_OHCU_RWDO_WNPRC','LPVGTI_OHCU_SKDS_CNPRC','LPVGTI_OHCU_ZJ_YNPRC','MK_TNPRC_JH_CPRC','MK_TNPRC_RWDO_WNPRC','MK_TNPRC_SKDS_CNPRC','RWNK_WNPRC_BF_ONPRC','RWNK_WNPRC_EV_NEPRC','RWNK_WNPRC_JH_CPRC','RWNK_WNPRC_LPVGTI_OHCU','RWNK_WNPRC_MK_TNPRC','RWNK_WNPRC_RWDO_WNPRC','RWNK_WNPRC_SKDS_CNPRC','RWNK_WNPRC_ZJ_YNPRC','SKDS_CNPRC_JH_CPRC','SKDS_CNPRC_RWDO_WNPRC','ZJ_YNPRC_BF_ONPRC','ZJ_YNPRC_JH_CPRC','ZJ_YNPRC_MK_TNPRC','ZJ_YNPRC_RWDO_WNPRC','ZJ_YNPRC_SKDS_CNPRC')
)


for (size in word) {
    for (i in seq(1,6)) {
        speciy <- species[i]
        pairses <- pairses_list[[i]]
        
        write_file_path <- file.path(root, "fst_out", paste(speciy, size,"fst_bglm.md", sep = "_"))
        write_file <- file(write_file_path, encoding = "UTF-8")
        sink(write_file,append = T)        

        for (pair in pairses) {

            filpath <- file.path( root,'fst_out',size,speciy,paste0("fst_",pair,'.density.out'))
            assign(paste(pair,'fst',sep='.'),read.table(filpath, header = T))
            datass <- eval(parse(text = paste(pair,'fst',sep='.')))[eval(parse(text = paste(pair,'fst',sep='.')))$containsgenes == 1,]
            datass <- na.exclude(datass)
            
            mod <- block.glm(c("Fst"),
                            c("containsNinteractgenez * number_of_genes"),
                            data=datass,
                            order.by = c("pos"),
                            strat.by = c("chr"),
                            block.size = 30, B = 1000)
            
            cat(paste("#",pair))
            print(as.data.frame(mod))
            cat("\n")
            
        }
        sink()
        
    }
}

# For Pi
pairses_list <- list(
    #c('aureus','fascicularis','thibetana_assamensis'),
    c('aureus','fascicularis','thibetana','assamensis'),
    c('Blue','Red','Purple'),
    c('M.brevicaudus','M.lasiotis','M.littoralis','M.mulatta','M.tcheliensis'),
    c("arctoides","assamensis","thibetana","mulata","fascicularis"),
    c('orange','brown','red'),
    c('JC_UoC','RWNK_WNPRC','LPVGTI_OHCU','EV_NEPRC','ZJ_YNPRC','BF_ONPRC','MK_TNPRC','SKDS_CNPRC','JH_CPRC','RWDO_WNPRC')
)

for (size in word) {

    for (i in seq(1,6)) {
        speciy <- species[i]
        pairses <- pairses_list[[i]]
        
        write_file_path <- file.path(root, "pi_out", paste(speciy, size,"pi_bglm.md", sep = "_"))
        write_file <- file(write_file_path, encoding = "UTF-8")
        sink(write_file,append = T)
        for (pair in pairses) {
            filpath <- file.path( root,'pi_out',size,speciy,paste0("pi_",pair,'.density.out'))
            assign(paste(pair,'pi',sep='.'),read.table(filpath, header = T))
            datass <- eval(parse(text = paste(pair,'pi',sep='.')))[eval(parse(text = paste(pair,'pi',sep='.')))$containsgenes == 1,]
            datass <- na.exclude(datass)
            
            mod <- block.glm(c("pi"),
                             c("containsNinteractgenez * number_of_genes"),
                             data=datass,
                             order.by = c("pos"),
                             strat.by = c("chr"),
                             block.size = 30, B = 1000)
            
            cat(paste("#",pair))
            print(as.data.frame(mod))
            cat("\n")
        }
        sink()
    }
}
