library(ggplot2)
library(pointr)
library(ggrepel)

options(scipen = 999)

setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/")
root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/'

type <- c("fst_out","pi_out")
word <- c("100k","30k")

########### For Fst

species <- c("chinese_mito","indian_mito","fas","aureus_new")
# ,,"chinese_org",,,"indian_org"
pairses_list <- list(
  c('Blue_Purple','Blue_Red','Purple_Red'),
  c('brown_green','brown_special','orange_special'),
  c ("arctoides_assamensis","arctoides_thibetana","arctoides_mulata","arctoides_fascicularis","assamensis_thibetana","assamensis_mulata","assamensis_fascicularis","thibetana_mulata","thibetana_fascicularis","mulata_fascicularis"),
  c('aureus_assamensis','aureus_thibetana','aureus_fascicularis','fascicularis_assamensis','fascicularis_thibetana','assamensis_thibetana')
  # c('M.brevicaudus_M.lasiotis','M.brevicaudus_M.mulatta','M.lasiotis_M.mulatta',
  #   'M.littoralis_M.brevicaudus','M.littoralis_M.lasiotis','M.littoralis_M.mulatta',
  #   'M.tcheliensis_M.brevicaudus','M.tcheliensis_M.lasiotis','M.tcheliensis_M.littoralis',
  #   'M.tcheliensis_M.mulatta'),
  # c('BF_ONPRC_JH_CPRC','BF_ONPRC_MK_TNPRC','BF_ONPRC_RWDO_WNPRC','BF_ONPRC_SKDS_CNPRC','EV_NEPRC_BF_ONPRC','EV_NEPRC_JH_CPRC','EV_NEPRC_MK_TNPRC','EV_NEPRC_RWDO_WNPRC','EV_NEPRC_SKDS_CNPRC','EV_NEPRC_ZJ_YNPRC','JC_UoC_BF_ONPRC','JC_UoC_EV_NEPRC','JC_UoC_JH_CPRC','JC_UoC_LPVGTI_OHCU','JC_UoC_MK_TNPRC','JC_UoC_RWDO_WNPRC','JC_UoC_RWNK_WNPRC','JC_UoC_SKDS_CNPRC','JC_UoC_ZJ_YNPRC','JH_CPRC_RWDO_WNPRC','LPVGTI_OHCU_BF_ONPRC','LPVGTI_OHCU_EV_NEPRC','LPVGTI_OHCU_JH_CPRC','LPVGTI_OHCU_MK_TNPRC','LPVGTI_OHCU_RWDO_WNPRC','LPVGTI_OHCU_SKDS_CNPRC','LPVGTI_OHCU_ZJ_YNPRC','MK_TNPRC_JH_CPRC','MK_TNPRC_RWDO_WNPRC','MK_TNPRC_SKDS_CNPRC','RWNK_WNPRC_BF_ONPRC','RWNK_WNPRC_EV_NEPRC','RWNK_WNPRC_JH_CPRC','RWNK_WNPRC_LPVGTI_OHCU','RWNK_WNPRC_MK_TNPRC','RWNK_WNPRC_RWDO_WNPRC','RWNK_WNPRC_SKDS_CNPRC','RWNK_WNPRC_ZJ_YNPRC','SKDS_CNPRC_JH_CPRC','SKDS_CNPRC_RWDO_WNPRC','ZJ_YNPRC_BF_ONPRC','ZJ_YNPRC_JH_CPRC','ZJ_YNPRC_MK_TNPRC','ZJ_YNPRC_RWDO_WNPRC','ZJ_YNPRC_SKDS_CNPRC')
)

size = "100k"

total_length <- 0
for (item in pairses_list) {
  total_length <- total_length + length(item)
}

est <- rep(0,total_length)
upper <- rep(0.0 ,total_length)
downer <- rep(0.0 ,total_length)
significance <- rep("" ,total_length)
total_lable <- rep("" ,total_length)

par(mfrow=c(5,4))  
cnt <- 0
for (i in seq(1,4)) {
  speciy <- species[i]
  pairses <- pairses_list[[i]]

    for (pairs in pairses) {
      cnt <- cnt + 1

      filpath <- file.path( root,'fst_out',size,speciy,paste0("fst_",pairs,'.density.out'))
      #assign(paste(pair,'fst',sep='.'),read.table(filpath, header = T))
      assign("pair",read.table(filpath, header = T))
      # ptr("pair",paste(pair,'fst',sep='.'))
      pair <- pair[pair$containsgenes == 1,]
      pair <- na.exclude(pair)
      
      
      mod <- lm(Fst ~ containsNinteractgenez * number_of_genes, data=pair)

      res1 <- mod

      est[cnt] <- res1$coefficients[2,1]
      significance[cnt] <- res1$coefficients[2,4]
      upper[cnt] <- confint(mod)[2,2]
      downer[cnt] <- confint(mod)[2,1]
      total_lable[cnt] <- pairs
      
      #       mod <- lm(Fst ~ containsNinteractgenez + number_of_genes, data=pair)
      #       acf(mod$residuals, main = pairs)
      # 
      #     }
      # }
        
      #   print(ggplot(pair, aes(x= seq(1:length(Fst)), y=sort(Fst)))+ geom_point() +
      #     geom_hline(yintercept = mean(pair$Fst), color = "red") +
      #     geom_text(aes(0,mean(pair$Fst),label = round(mean(pair$Fst), digits = 4), vjust = -1), color = 'red') +
      #     xlab(NULL) + ylab("Fst") + ggtitle(pairs)
      #   )
      # }
      # plot(y=sort(pair$Fst), x = seq(1:length(pair$Fst)),main=pairs)
        #hist(pair$Fst, main = pairs, breaks = seq(-0.1,0.2,0.005))
    }
}

  values <- round(est, digits = 4)
  for (index in seq(1,length(significance))) {
    if (significance[index] < 0.05 ) {
      if ( 0 < downer[index] | 0 > upper[index]  ) {
        values[index] <- paste0(values[index],"*")
      }
    }
  }
  
  for_plot <- data.frame(
    est,
    upper,
    downer,
    values
  ) 
  
  total_lable[11] <- "assamensis_thibetana_1"
  
  for_plot$lable <- factor(
    total_lable,
    rev(total_lable),
    rev(c("Blue:Purple","Blue:Red","Purple:Red","Brown:Orange","Brown:Special","Orange:Special","M.arctoides:M.assamensis",
          "M. arctoides:M. thibetana","M. arctoides:M. mulatta","M. arctoides:M. fascicularis","M. assamensis_1:M. thibetana",
          "M. assamensis:M. mulatta","M. assamensis:M. fascicularis","M. thibetana:M. mulatta","M. thibetana:M. fascicularis",
          "M. mulatta:M. fascicularis",
          "M. f. aurea:M. assamensis","M. f. aurea:M. thibetana","M. f. aurea:M. fascicularis", 'M. fascicularis:M. assamensis','M. fascicularis:M. thibetana', 'M. assamensis_2:M. thibetana'))
  )
  c('aureus_assamensis','aureus_thibetana','aureus_fascicularis','fascicularis_assamensis','fascicularis_thibetana','assamensis_thibetana')
  
  interested <- c("Blue:Purple","Blue:Red","Purple:Red","Brown:Orange","Brown:Special","Orange:Special",
                  "M. arctoides:M. thibetana","M. arctoides:M. mulatta","M. arctoides:M. fascicularis","M.arctoides:M.assamensis",
                  "M. f. aurea:M. fascicularis","M. f. aurea:M. thibetana","M. f. aurea:M. assamensis")
  a <- rev(ifelse(for_plot$lable %in% interested, "red", "black"))
  
  
  pdf(paste0("fst_",size,"_onlygenic.pdf"))
  ggplot(for_plot, aes(x = lable, y= round(est, digits = 4))) + 
    geom_hline(yintercept = 0, color = "red") + 
    geom_pointrange(aes(ymax = upper, ymin=downer), color = "darkblue") + 
    geom_text(aes(label=values),nudge_x=0.3,nudge_y=-0.005) +
    scale_x_discrete("", labels = rev(
      c("Blue : Purple",expression("Blue"~":"~"Red"["China"]),expression("Purple"~":"~"Red"["China"]),
        "Brown : Orange",expression("Brown"~":"~"Red"["India"]),expression("Orange"~":"~"Red"["India"]),
        expression(italic("M.arctoides")~":"~italic(" M.assamensis"[1])),expression(italic("M. arctoides")~":"~italic("M. thibetana"[1])),
        expression(italic("M. arctoides")~":"~italic("M. mulatta")),expression(italic("M. arctoides")~":"~italic("M. fascicularis")[1]),expression(italic("M. assamensis"[1])~":"~italic("M. thibetana"[1])),
        expression(italic("M. assamensis"[1])~":"~italic("M. mulatta")),expression(italic("M. assamensis"[1])~":"~italic("M. fascicularis"[1])),expression(italic("M. thibetana"[1])~":"~italic("M. mulatta")),expression(italic("M. thibetana"[1])~":"~italic("M. fascicularis"[1])),
        expression(italic("M. mulatta")~":"~italic("M. fascicularis"[1])),
        expression(italic("M. f. aurea")~":"~italic("M. assamensis"[2])),
        expression(italic("M. f. aurea")~":"~italic("M. thibetana"[2])),
        expression(italic("M. f. aurea")~":"~italic("M. fascicularis"[2])),
        expression(italic("M. fascicularis"[2])~":"~italic("M. assamensis"[2])),
        expression(italic("M. fascicularis"[2])~":"~italic("M. thibetana"[2])),
        expression(italic("M. assamensis"[2])~":"~italic("M. thibetana"[2]))
      )
    )) + 
    #scale_y_continuous(limits = c(-0.025, 0.075), breaks = seq(-0.025,0.075, 0.025)) + 
    theme_bw() + 
    theme(text = element_text(size = 10),axis.text=element_text(size=12),axis.text.y = element_text(hjust=1, face = "italic",colour = a)) + 
    ylab(NULL) + 
    coord_flip()
  dev.off()
  
############# For Pi

species <- c("chinese_mito","indian_mito","fas","aureus_new")
pairses_list <- list(
  c('Blue','Red','Purple'),
  c('orange','brown','red'),
  c("arctoides","assamensis","thibetana","mulata","fascicularis"),
  c('aureus','fascicularis','thibetana','assamensis')
)

total_length <- 0
for (item in pairses_list) {
  total_length <- total_length + length(item)
}
est <- rep(0,total_length)
upper <- rep(0.0 ,total_length)
downer <- rep(0.0 ,total_length)
significance <- rep("" ,total_length)
total_lable <- rep("" ,total_length)


size = "100k"
cnt <- 0
for (i in seq(1,4)) {
  speciy <- species[i]
  pairses <- pairses_list[[i]]
  
  for (pair in pairses) {
    cnt <- cnt + 1
    
    filpath <- file.path( root,'pi_out',size,speciy,paste0("pi_",pair,'.density.out'))
    assign("datas",read.table(filpath, header = T))
    # ptr(pair,paste(pair,'pi',sep='.'))
    datas <- datas[datas$containsgenes == 1,]
    mod <- lm(pi ~ containsNinteractgenez * number_of_genes, data=datas)
    
    
    res1 <- summary(mod)
    
    est[cnt] <- res1$coefficients[2,1]
    significance[cnt] <- res1$coefficients[2,4]
    upper[cnt] <- confint(mod)[2,2]
    downer[cnt] <- confint(mod)[2,1]
    total_lable[cnt] <- pair
  }
}

values <- round(est, digits = 4)
for (index in seq(1,length(significance))) {
  if (significance[index] < 0.05 ) {
    if ( 0 < downer[index] | 0 > upper[index]  ) {
      values[index] <- paste0(values[index],"*")
    }
  }
}

for_plot <- data.frame(
  est,
  upper,
  downer,
  values
) 

total_lable[8] <-'assamensis_1'
total_lable[9] <- 'thibetana_2'
total_lable[11] <- "fascicularis_1"
total_lable[13] <- "fascicularis_2"

for_plot$lable <- factor(
  total_lable,
  rev(total_lable),
  rev(c("Blue","Red","Purple","Orange","Brown",'Special',"M. arctoides","M. assamensis_1", "M. thibetana_1", 
        "M. mulatta","M. fascicularis_1", "M. f. aurea", "M. fascicularis_2", "M. thibetana","M. assamensis"))
)

interested <- c("Red", "Blue", "Green", 'Special',"Orange","Purple", "Brown", "M. arctoides", "M. f. aurea")
a <- rev(ifelse(for_plot$lable %in% interested, "red", "black"))

pdf(paste0("pi_",size,"_onlygenic.pdf"))
ggplot(for_plot, aes(x = lable, y= round(est, digits = 4))) + 
  geom_hline(yintercept = 0, color = "red") + 
  geom_pointrange(aes(ymax = upper, ymin=downer), color = "darkblue") + 
  geom_text(aes(label=values),nudge_x=0.3,nudge_y=-0.005) +
  scale_x_discrete("", labels=  rev(c("Blue",expression("Red"["China"]),"Purple","Orange","Brown",expression('Red'['India']),expression(italic("M. arctoides")),expression(italic("M. assamensis"[1])), expression(italic("M. thibetana"[1])), 
                                      expression(italic("M. mulatta")),expression(italic("M. fascicularis"[italic(1)])), expression(italic("M. f. aurea")), expression(italic("M. fascicularis"[2])), expression(italic("M. thibetana"[2])),expression(italic("M. assamensis"[2]))  ))) + 
  # scale_y_continuous(limits = c(-0.025, 0.01), breaks = seq(-0.025,0.01, 0.01)) + 
  theme_bw() + 
  theme(text = element_text(size = 10),axis.text=element_text(size=12),axis.text.y = element_text(hjust=1, face = "italic", colour = a)) + 
  ylab(NULL) + 
  coord_flip()
dev.off()

############# For ROH

species <- c("chinese_mito","indian_mito","fas","aureus_new")
pairses_list <- list(
  c('Blue','Red','Purple'),
  c('orange','brown','red'),
  c("arctoides","assamensis","thibetana","mulata","fascicularis"),
  c('aureus','fascicularis','thibetana','assamensis')
)

total_length <- 0
for (item in pairses_list) {
  total_length <- total_length + length(item)
}
est <- rep(0,total_length)
upper <- rep(0.0 ,total_length)
downer <- rep(0.0 ,total_length)
significance <- rep("" ,total_length)
total_lable <- rep("" ,total_length)


cnt <- 0
for (i in seq(1,4)) {
  speciy <- species[i]
  pairses <- pairses_list[[i]]
  
  for (pair in pairses) {
    cnt <- cnt + 1
    
    filpath <- file.path( root,'roh_out',speciy,paste0("roh_",pair,'.density.out'))
    assign(paste(pair,'roh',sep='.'),read.table(filpath, header = T))
    # ptr(pair,paste(pair,'roh',sep='.'))
    ptr("now_data", paste(pair,'roh',sep='.'))
    
    now_data <- now_data[(now_data$containsgenes == 1),]
    now_data$length <- now_data$length/100000
    now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 
    
    mod <- lm(log(length) ~ containsNinteractgenez*num_genes, data=now_data)
    
    res1 <- summary(mod)
    confint(mod)
    
    est[cnt] <- res1$coefficients[2,1]
    significance[cnt] <- res1$coefficients[2,4]
    upper[cnt] <- confint(mod)[2,2]
    downer[cnt] <- confint(mod)[2,1]
    total_lable[cnt] <- pair
  }
}

setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/roh_results")

pairses_list <- list(
    c('Blue','Red','Purple'),
    c('orange','brown','red'),
    c("arctoides","assamensis","thibetana","mulata","fascicularis"),
    c('aureus','fascicularis','thibetana','assamensis')
)


for (i in seq(1,4)) {
    speciy <- species[i]
    pairses <- pairses_list[[i]]
    write_file <- file(paste0("lm_roh_",speciy,".md"), encoding = "UTF-8")
    sink(write_file,append = T)
    for (pair in pairses) {
        
        filpath <- file.path( root,'roh_out',speciy,paste0("roh_",pair,'.density.out'))
        assign(paste(pair,'roh',sep='.'),read.table(filpath, header = T))
        # ptr(pair,paste(pair,'roh',sep='.'))
        ptr("now_data", paste(pair,'roh',sep='.'))
        
        now_data <- now_data[(now_data$containsgenes == 1),]
        now_data$length <- now_data$length/100000
        now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 
        
        mod <- lm(log(length) ~ containsNinteractgenez*num_genes, data=now_data)
        
        # est <- mod$coefficients
        # CI <- confint(mod)[,1]
        ci = coef(mod) - confint(mod)[,1]
        
        cat(paste("##",pair))
        print(data.frame("est" = coef(mod), "conf" = ci))
        
    }
    sink()
}


values <- round(est, digits = 4)
for (index in seq(1,length(significance))) {
  if (as.numeric(significance[index]) < 0.05 ) {
    if ( 0 < downer[index] | 0 > upper[index]  ) {
      values[index] <- paste0(values[index],"*")
    }
  }
}

for_plot <- data.frame(
  est,
  upper,
  downer,
  values
) 

total_lable[8] <-'assamensis_1'
total_lable[9] <- 'thibetana_2'
total_lable[11] <- "fascicularis_1"
total_lable[13] <- "fascicularis_2"

for_plot$lable <- factor(
  total_lable,
  rev(total_lable),
  rev(c("Blue","Red","Purple","Orange","Brown",'Special',"M. arctoides","M. assamensis_1", "M. thibetana_1", 
        "M. mulatta","M. fascicularis_1", "M. f. aurea", "M. fascicularis_2", "M. thibetana","M. assamensis"))
)

interested <- c("Red", "Blue", "Green", "Orange","Purple", "Special","Brown", "M. arctoides", "M. f. aurea")
a <- rev(ifelse(for_plot$lable %in% interested, "red", "black"))

pdf("roh.pdf")
ggplot(for_plot, aes(x = lable, y= round(est, digits = 4))) + 
  geom_hline(yintercept = 0, color = "red") + 
  geom_pointrange(aes(ymax = upper, ymin=downer), color = "darkblue") + 
  geom_text(aes(label=values),nudge_x=0.3,nudge_y=-0.005) +
  scale_x_discrete("", labels=  rev(c("Blue",expression("Red"["China"]),"Purple","Orange","Brown",expression('Red'['India']),expression(italic("M. arctoides")),expression(italic("M. assamensis"[1])), expression(italic("M. thibetana"[1])), 
                                      expression(italic("M. mulatta")),expression(italic("M. fascicularis"[italic(1)])), expression(italic("M. f. aurea")), expression(italic("M. fascicularis"[2])), expression(italic("M. thibetana"[2])),expression(italic("M. assamensis"[2]))  ))) + 
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0,1.5,0.5)) + 

  theme_bw() + 
  theme(text = element_text(size = 10),axis.text=element_text(size=12),axis.text.y = element_text(hjust=1, colour = a)) + 
  ylab(NULL) + 
  coord_flip()
dev.off()


