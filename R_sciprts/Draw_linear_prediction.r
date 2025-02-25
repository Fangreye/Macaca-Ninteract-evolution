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
library(ggiraphExtra)
library(ggggeffects)
library(stringr)

show_percent <- function(x) {  
    perc <- sprintf("%0.01f%%",x*100) 
    return(perc)
}
scaleFUN <- function(x) sprintf("%.1f", x)

draw_this <- function(lmodel, x_limit, species_lable){
    # lmodel <- Blue_mod
    # x_limit = 5
    
    to_draw <- paste0("num_genes [1:", x_limit,"]")
    temp <- ggpredict(lmodel,  c(to_draw,"containsNinteractgenez"), back.transform = F)
    
    num_gene_var <- lmodel$model$num_genes
    num_gene_table <- as.data.frame(table(num_gene_var)/length(num_gene_var))
    num_gene_table <- num_gene_table[as.numeric(num_gene_table$num_gene_var) <= x_limit,]
    num_gene_table$Percent <- sapply(num_gene_table$Freq,FUN=show_percent)
    
    if (x_limit == 2) {
        pos_seq <- c(1.15, x_limit - 0.15)
    } else if (x_limit == 3) {
        pos_seq <- c(1.25, seq(2,x_limit-1),x_limit - 0.25)
    } else {
        pos_seq <- c(1.25, 2.10, seq(3,x_limit-1), x_limit - 0.25)
    }
    
    # if (x_limit <= 4) {
    #     lable_size = 8
    # } else if (x_limit <= 6) {
    #     lable_size = 7
    # } else {
    #     lable_size = 5
    # }
    # lable_size = 8 / ceiling(x_limit/2) 
    lable_size = 7
    #lable_pos = min(temp$predicted) + 0.01* abs(min(temp$predicted))
    lable_measure = max(temp$predicted) - min(temp$predicted)
    lable_pos = min(temp$predicted) + 0.15 * lable_measure
    
    something <- as.data.frame(temp)
    ggplot(something, aes(x = x, y = predicted, colour = group)) + 
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), alpha = .15) + 
        annotate("text", x=pos_seq, y=lable_pos, label= num_gene_table$Percent, size = lable_size, angle = 45) +
        theme(panel.background = element_blank(), legend.position = "none",
              panel.border = element_rect(fill = NA,  colour = "grey20"),
              axis.title.x = element_text(size=22), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ,
              #plot.margin = unit(rep(10,4), "pt"), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        ) + 
        scale_colour_grey(start = 0.2, end = 0.6) + 
        #scale_color_manual(values=c('blue','red')) + 
        #scale_color_manual(values=c('gray30','gray10')) +
        scale_fill_grey(start = 0.2, end = 0.6) + 
        #scale_linetype_manual(values=c("dash", "dotted")) + 
        # scale_fill_manual(values=c('blue', 'red'), name="fill") +
        # scale_fill_manual(values=c('gray30', 'gray10'), name="fill") +
        scale_x_continuous(name=species_lable, limits=c(1,x_limit), breaks = 1:x_limit) +
        scale_y_continuous(labels=scaleFUN) +
        labs(x = species_lable, y = "")
}


setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/roh_results/")

theme_set(theme_ggeffects())

# species_names <- c("aureus","fascicularis","thibetana","assamensis")
# species_names <- c("arctoides","assamensis","thibetana","mulata","fascicularis")
# species_names <- c('Blue','Red','Purple')
# species_names <- c('orange','brown','red')

species_names <- list(
    c("aureus","fascicularis","thibetana","assamensis"),
    c("arctoides","assamensis","thibetana","mulata","fascicularis"),
    c('Blue','Red','Purple'),
    c('orange','brown','red')
)
place <- c(
    "aureus_new","fas","chinese_mito","indian_mito"
)
# root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/aureus_new/"
# root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/fas/"
# root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/chinese_mito/"
# root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/indian_mito/"

root <- "F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/roh_out/"

for (index in 1:length(place)){
    
    dataset <- place[index]
    species_data <- species_names[[index]]
    
    for (species in species_data) {
        assign(paste0(dataset,species),read.table(file.path(root,dataset,paste0("roh_",species,".density.out")), header = T)) 
        ptr("now_data", paste0(dataset,species))
        now_data <- now_data[now_data$containsgenes == 1,]
        now_data$length <- now_data$length/100000
        now_data$log_length <- log(now_data$length)
        now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 
        
        mean(now_data$length[now_data$num_genes == 2])
        assign(paste0(dataset,species,"_mod"),lm(log_length ~ containsNinteractgenez * num_genes, data= now_data))
    }
    
    # read.table(file.path(root,'ROH_Blue_roh_density.out'), header=T)
    # assign(species,read.table(file.path(root,paste0("roh_",species,".density.out")), header = T)) 
    # ptr("now_data", species)
    # now_data <- now_data[now_data$containsgenes == 1,]
    # now_data$length <- now_data$length/100000
    # now_data$log_length <- log(now_data$length)
    # now_data$gene_density_on_ROHs <- (now_data$num_genes/now_data$length) 
    # 
    # mean(now_data$length[now_data$num_genes == 2])
    # 
    # assign(paste0(species,"_mod"),lm(log_length ~ containsNinteractgenez * num_genes, data= now_data))
}

species_names <- c('Blue','Red','Purple')
species_names <- c('orange','brown','red')
species_names <- c("arctoides","assamensis","thibetana","mulata","fascicularis")
species_names <- c("aureus","fascicularis","thibetana","assamensis")
labels <- list(
    c("f. aurea","fascicularis","thibetana","assamensis"),
    c("arctoides","assamensis","thibetana","mulatta","fascicularis"),
    c('Blue','Red','Purple'),
    c('orange','brown','red')
)
for (index in 1:length(place)) {
    
    dataset <- place[index]
    species_data <- species_names[[index]]
    label_here <- labels[[index]]
    
    if (index %in% 1:2) {
        for (i in 1:length(species_data)) {
            species = species_data[i]
            
            my_label = label_here[i]
            a = eval(parse(text=paste0(dataset,species,"_mod")))
            x_limit = quantile(a$model$num_genes, probs = seq(0.9, 0.9))
            assign(paste0(dataset,species,"_plot"), draw_this(eval(parse(text=paste0(dataset,species,"_mod"))),x_limit, 
                                                            as.expression(bquote(italic(.(paste0("M. ",my_label)))))
                                                            
                                                            #lmodel = aureus_newassamensis_mod 
                                                            # species_lable = as.expression(bquote(italic(.(paste0("M. ",my_label)))))
            ))
        }
    } else {
        for (i in 1:length(species_data)) {
            species = species_data[i]
            my_label = label_here[i]
            
            a = eval(parse(text=paste0(dataset,species,"_mod")))
            x_limit = quantile(a$model$num_genes, probs = seq(0.9, 0.9))
            assign(paste0(dataset,species,"_plot"), draw_this(eval(parse(text=paste0(dataset,species,"_mod"))),x_limit, 
                                                            str_to_title(my_label))
            )
        }
        
    }
}
        

bottom_text = text_grob("Number of genes",size = 20,just = 'center')
# bottom_text = text_grob("gene_density",size = 20,just = 'center')
bottom_text$vp <- viewport(x=0.5, y=1)

left_text = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 1,size = 20)
left_text$vp <- viewport(x=0.2, y=0.2)

chinese_grid <- arrangeGrob(
    grobs = list(chinese_mitoBlue_plot,chinese_mitoRed_plot,chinese_mitoPurple_plot),
    ncol = 4
)

indian_grid <- arrangeGrob(
    grobs = list(indian_mitoorange_plot,indian_mitobrown_plot,indian_mitored_plot),
    ncol = 4
)

arctoides_grid <- arrangeGrob(
    grobs = list(fasarctoides_plot,fasassamensis_plot,fasthibetana_plot,fasmulata_plot,fasfascicularis_plot),
    ncol = 4
)

aureus_grid <- arrangeGrob(
    grobs = list(aureus_newaureus_plot,aureus_newfascicularis_plot,aureus_newthibetana_plot,aureus_newassamensis_plot),
    ncol = 4
)

something = rbind(chinese_grid, indian_grid, arctoides_grid, aureus_grid)

setwd("F:/Work/Now_work/2021_05_18.ben_/06.corrected/predicts")
png("Predicted_values_ROH_gray.png", width = 480, height = 360, units='mm', res = 300)   
grid.arrange(
    #arrangeGrob(grobs = list(chinese_grid,indian_grid,arctoides_grid), ncol = 1),
    something,
    left = text_grob("ln(ROH length (100 kb))", rot = 90, vjust = 0.5,hjust = 0.5,size = 20),
    bottom = bottom_text
)
dev.off()

###################
scaleFUN <- function(x) sprintf("%.3f", x)

draw_this <- function(lmodel, x_limit, species_lable){
    # lmodel <- Blue_mod
    # x_limit = 5
    
    to_draw <- paste0("number_of_genes [1:", x_limit,"]")
    temp <- ggpredict(lmodel,  c(to_draw,"containsNinteractgenez"), back.transform = F)
    
    num_gene_var <- lmodel$model$number_of_genes
    num_gene_table <- as.data.frame(table(num_gene_var)/length(num_gene_var))
    num_gene_table <- num_gene_table[as.numeric(num_gene_table$num_gene_var) <= x_limit,]
    num_gene_table$Percent <- sapply(num_gene_table$Freq,FUN=show_percent)
    
    if (x_limit == 2) {
        pos_seq <- c(1.15, x_limit - 0.15)
    } else if (x_limit == 3) {
        pos_seq <- c(1.25, seq(2,x_limit-1),x_limit - 0.25)
    } else {
        pos_seq <- c(1.25, 2.10, seq(3,x_limit-1), x_limit - 0.25)
    }
    
    # if (x_limit <= 4) {
    #     lable_size = 8
    # } else if (x_limit <= 6) {
    #     lable_size = 7
    # } else {
    #     lable_size = 5
    # }
    # lable_size = 8 / ceiling(x_limit/2) 
    lable_size = 7
    #lable_pos = min(temp$predicted) + 0.01* abs(min(temp$predicted))
    lable_measure = max(temp$conf.high) - min(temp$conf.low)
    # if (lable_measure < )
    
    lable_pos = min(temp$conf.low) + 0.2 * lable_measure
    # print(paste(lable_measure,":", min(temp$conf.low)))
    
    something <- as.data.frame(temp)
    ggplot(something, aes(x = x, y = predicted, colour = group)) + 
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), alpha = .15) + 
        annotate("text", x=pos_seq, y=lable_pos, label= num_gene_table$Percent, size = lable_size, angle = 45) +
        theme(panel.background = element_blank(), legend.position = "none",panel.border = element_rect(fill = NA,  colour = "grey20"),
              axis.title.x = element_text(size=22), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ,
              #plot.margin = unit(rep(10,4), "pt"), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        ) + 
        scale_colour_grey(start = 0.2, end = 0.6) + 
        scale_fill_grey(start = 0.2, end = 0.6) + 
        # scale_color_manual(values=c('blue','red')) + 
        # scale_fill_manual(values=c('blue', 'red'), name="fill") +
        scale_x_continuous(name=species_lable, limits=c(1,x_limit), breaks = 1:x_limit) +
        scale_y_continuous(labels=scaleFUN) +
        labs(x = species_lable, y = "")
}


species_names <- list(
    c('Blue','Red','Purple'),
    c('orange','brown','red'),
    c("arctoides","assamensis","thibetana","mulata","fascicularis"),
    c("aureus","fascicularis","thibetana","assamensis")
)
place <- c(
    "chinese_mito","indian_mito","fas","aureus_new"
)

lables <- list(
    c("Blue","Red","Purple"),
    c("Orange","Brown",'Red') ,
    c(
        as.expression(bquote(italic("M. arctoides"))), as.expression(bquote(italic("M. assamensis"))),
        as.expression(bquote(italic("M. thibetana"))), as.expression(bquote(italic("M. mulatta"))),
        as.expression(bquote(italic("M. fascicularis")))
    ),
    c(
        as.expression(bquote(italic("M. f. aurea"))), as.expression(bquote(italic("M. fascicularis"))), 
        as.expression(bquote(italic("M. thibetana"))),as.expression(bquote(italic("M. assamensis")))
    )
)

as.expression(bquote(italic("A")))

root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/'

size = "30k"
for (index in seq(1,4)) {

    dataset <- place[index]
    species_data <- species_names[[index]]
    this_lable <- lables[[index]]
    

    sometext <- rep("",length(species_data))
    
    for (sp_i in 1:length(species_data)) {
        speciy <- species_data[sp_i]
        labmy_labelle_here <- this_lable[sp_i]
        
        this_file <- file.path( root,'pi_out',size,dataset,paste0("pi_",speciy,'.density.out'))
        
        assign(paste0(dataset,speciy),read.table(this_file, header = T)) 
        ptr("now_data", paste0(dataset,speciy))
        now_data <- now_data[now_data$containsgenes == 1,]
        
        assign(paste0(dataset,speciy,"_mod"),lm(pi ~ containsNinteractgenez * number_of_genes, data= now_data))
        
        x_limit = quantile(now_data$number_of_genes, probs = seq(0.9, 0.9))
        assign(paste0(dataset,speciy,"_plot"), draw_this(eval(parse(text=paste0(dataset,speciy,"_mod"))),x_limit, 
                                                          labmy_labelle_here
                                                          )
        )
        
        sometext[sp_i] <- paste0(dataset,speciy,"_plot")
    }
    cmd <- paste0("arrangeGrob(
        grobs = list(",paste(sometext, collapse = ","),"),
        ncol = 4
    )")
    
    assign(paste0(dataset,"_grob"), eval(parse(text = cmd)))
}

something = rbind(chinese_mito_grob, indian_mito_grob, fas_grob, aureus_new_grob)
png(paste0("Predicted_values_","pi","_",size,"_gray.png"), width = 480, height = 360, units='mm', res = 300)   
grid.arrange(
    #arrangeGrob(grobs = list(chinese_grid,indian_grid,arctoides_grid), ncol = 1),
    something,
    left = text_grob(expression(pi), rot = 90, vjust = 0.5,hjust = 0.5,size = 30),
    bottom = bottom_text
)
dev.off()

###################

scaleFUN <- function(x) sprintf("%.3f", x)

draw_this <- function(lmodel, x_limit, species_lable){
    
    to_draw <- paste0("number_of_genes [1:", x_limit,"]")
    temp <- ggpredict(lmodel,  c(to_draw,"containsNinteractgenez"), back.transform = F)
    
    num_gene_var <- lmodel$model$number_of_genes
    num_gene_table <- as.data.frame(table(num_gene_var)/length(num_gene_var))
    num_gene_table <- num_gene_table[as.numeric(num_gene_table$num_gene_var) <= x_limit,]
    num_gene_table$Percent <- sapply(num_gene_table$Freq,FUN=show_percent)
    
    if (x_limit == 2) {
        pos_seq <- c(1.15, x_limit - 0.15)
    } else if (x_limit == 3) {
        pos_seq <- c(1.25, seq(2,x_limit-1),x_limit - 0.25)
    } else {
        pos_seq <- c(1.25, 2.10, seq(3,x_limit-1), x_limit - 0.25)
    }
    
    lable_size = 7
    lable_measure = max(temp$conf.high) - min(temp$conf.low)
    lable_pos = min(temp$conf.low) + 0.2 * lable_measure
    # print(paste(lable_measure,":", min(temp$conf.low)))
    
    something <- as.data.frame(temp)
    ggplot(something, aes(x = x, y = predicted, colour = group)) + 
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), alpha = .15) + 
        annotate("text", x=pos_seq, y=lable_pos, label= num_gene_table$Percent, size = lable_size, angle = 45) +
        theme(panel.background = element_blank(), legend.position = "none",panel.border = element_rect(fill = NA,  colour = "grey20"),
              axis.title.x = element_text(size=21), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ,
              #plot.margin = unit(rep(10,4), "pt"), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        ) + 
        scale_colour_grey(start = 0.2, end = 0.6) + 
        scale_fill_grey(start = 0.2, end = 0.6) + 
        #scale_color_manual(values=c('blue','red')) + 
        #scale_fill_manual(values=c('blue', 'red'), name="fill") +
        scale_x_continuous(name=species_lable, limits=c(1,x_limit), breaks = 1:x_limit) +
        scale_y_continuous(labels=scaleFUN) +
        labs(x = species_lable, y = "")
}


species_names <- list(
    c('Blue_Purple','Blue_Red','Purple_Red'),
    c('brown_orange','brown_red','orange_red'),
    c ("arctoides_assamensis","arctoides_thibetana","arctoides_mulata","arctoides_fascicularis","assamensis_thibetana","assamensis_mulata","assamensis_fascicularis","thibetana_mulata","thibetana_fascicularis","mulata_fascicularis"),
    c('aureus_assamensis','aureus_thibetana','aureus_fascicularis','fascicularis_assamensis','fascicularis_thibetana','assamensis_thibetana')
)
place <- c(
    "chinese_mito","indian_mito","fas","aureus_new"
)

lables <- list(
    c("Blue : Purple","Blue : Red","Purple : Red"),
    c("Brown : Orange","Brown : Red","Orange : Red" ),
    c(
        as.expression(bquote(italic("M.arctoides")~":"~italic(" M.assamensis"))),as.expression(bquote(italic("M. arctoides")~":"~italic("M. thibetana"))),as.expression(bquote(italic("M. arctoides")~":"~italic("M. mulatta"))),
        as.expression(bquote(italic("M. arctoides")~":"~italic("M. fascicularis"))),as.expression(bquote(italic("M. assamensis")~":"~italic("M. thibetana"))),as.expression(bquote(italic("M. assamensis")~":"~italic("M. mulatta"))),
        as.expression(bquote(italic("M. assamensis")~":"~italic("M. fascicularis"))),as.expression(bquote(italic("M. thibetana")~":"~italic("M. mulatta"))),as.expression(bquote(italic("M. thibetana")~":"~italic("M. fascicularis"))),
        as.expression(bquote(italic("M. mulatta")~":"~italic("M. fascicularis")))
    ),
    c(
        as.expression(bquote(italic("M. f. aurea")~":"~italic("M. assamensis"))),
        as.expression(bquote(italic("M. f. aurea")~":"~italic("M. thibetana"))),
        as.expression(bquote(italic("M. f. aurea")~":"~italic("M. fascicularis"))),
        as.expression(bquote(italic("M. fascicularis")~":"~italic("M. assamensis"))),
        as.expression(bquote(italic("M. fascicularis")~":"~italic("M. thibetana"))),
        as.expression(bquote(italic("M. assamensis")~":"~italic("M. thibetana")))
    )
)

root <- 'F:/Work/Now_work/2021_05_18.ben_/06.corrected/00.output/'

size = "30k"
for (index in seq(1,4)) {
    
    dataset <- place[index]
    species_data <- species_names[[index]]
    this_lable <- lables[[index]]
    
    
    sometext <- rep("",length(species_data))
    
    for (sp_i in 1:length(species_data)) {
        speciy <- species_data[sp_i]
        labmy_labelle_here <- this_lable[sp_i]
        
        this_file <- file.path( root,'fst_out',size,dataset,paste0("fst_",speciy,'.density.out'))
        
        assign(paste0(dataset,speciy),read.table(this_file, header = T)) 
        ptr("now_data", paste0(dataset,speciy))
        now_data <- now_data[now_data$containsgenes == 1,]
        
        assign(paste0(dataset,speciy,"_mod"),lm(Fst ~ containsNinteractgenez * number_of_genes, data= now_data))
        
        x_limit = quantile(now_data$number_of_genes, probs = seq(0.9, 0.9))
        assign(paste0(dataset,speciy,"_plot"), draw_this(eval(parse(text=paste0(dataset,speciy,"_mod"))),x_limit, 
                                                         labmy_labelle_here
        )
        )
        
        sometext[sp_i] <- paste0(dataset,speciy,"_plot")
    }
    cmd <- paste0("arrangeGrob(
        grobs = list(",paste(sometext, collapse = ","),"),
        ncol = 4)"
        #ncol = 4 , padding = unit(0.5, 'line'))")
    )
    
    assign(paste0(dataset,"_grob"), eval(parse(text = cmd)))
}

something = rbind(chinese_mito_grob, indian_mito_grob, fas_grob, aureus_new_grob)
png(paste0("Predicted_values_","fst","_",size,".png"), width = 540, height = 360/5*7, units='mm', res = 300)   
grid.arrange(
    #arrangeGrob(grobs = list(chinese_grid,indian_grid,arctoides_grid), ncol = 1),
    something,
    left = text_grob(expression(F[ST]), rot = 90, vjust = 0.5,hjust = 0.5,size = 30),
    bottom = bottom_text,
    right = text_grob("", rot = 90, vjust = 0.5,hjust = 0.5,size = 30),
    top = text_grob("",size = 20,just = 'center')
)
dev.off()


