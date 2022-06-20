N_interact <- read.csv("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/10.a_output/Chinese_Ninteract.csv", header = F)
Normal <- read.csv("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/10.a_output/Chinese_Normal.csv", header = F)
colnames(N_interact) <- c("gene","dN","dS")
colnames(Normal) <- c("gene","dN","dS")

# N_interact <- N_interact[N_interact$V3 != 0,]
# Normal <- Normal[Normal$V3 != 0,]

dNdS.Ninteract <- N_interact$dN / N_interact$dS
N_interact.dNdS <- data.frame(dNdS = dNdS.Ninteract, type = rep(1,length(N_interact[,1])))
N_interact <- cbind(N_interact, N_interact.dNdS)

N_interact.meaningful <- N_interact[!is.na(N_interact.dNdS$dNdS) & !is.infinite(N_interact.dNdS$dNdS),]
length(N_interact.meaningful$dNdS[N_interact.meaningful$dNdS<5])
# N_interact.dNdS <- N_interact.dNdS[N_interact.dNdS$dNdS < 5,]
dNdS.Normal <- Normal$dN / Normal$dS
Normal.dNdS <- data.frame(dNdS = dNdS.Normal, type = rep(0,length(Normal[,1])))
Normal <- cbind(Normal, Normal.dNdS)

Normal.meaningful <- Normal[!is.na(Normal.dNdS$dNdS) & !is.infinite(Normal.dNdS$dNdS),]
length(Normal.meaningful$dNdS[Normal.meaningful$dNdS<5])
# Normal.dNdS <- Normal.dNdS[Normal.dNdS$dNdS < 5,]

data1 <- rbind(N_interact,Normal)
data <- rbind(N_interact.meaningful, Normal.meaningful)

png("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/10.a_output/Chinese_stats.jpg", width = 1500, height = 900)
par(mfrow=c(1,4), cex.axis=2  ,cex.lab  = 3, mar = c(8, 8, 4, 4))
hist(data$dNdS[data$dNdS < 5], breaks = seq(0,5,0.25), main= paste0("Amount of samples:", length(data$dNdS[data$dNdS < 5]) , "/",length(data1$dNdS)))
plot(data1$dS,data1$dN, pch = 1, main=paste0('Meaningful dNdS: ', length(data$dNdS),'/',length(data1$dN)))
# with(data1[(data1$dN < 0.5) & (data1$dS < 0.5),], plot(dS, dN))
with(data[data$type ==0 ,], plot(dS, dN, main = paste0('Meaningful Normal:' , length(data$dN[data$type ==0]),'/', length(data1$dN[data1$type ==0]) )))
with(data[data$type ==1 ,], plot(dS, dN, main = paste0('Meaningful N-interact:' , length(data$dN[data$type ==1]) ,'/', length(data1$dN[data1$type ==1])) ))
dev.off()

# png("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/05.values/Chinese_cate.jpg", width = 1200, height = 700)
# par(mfrow=c(1,2), cex.axis=2  ,cex.lab  = 3, mar = c(8, 8, 4, 4))
# 
# dev.off()

# png("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/05.values/Chinese_stats.jpg", width = 1200, height = 800)
# par(mfrow=c(1,3))
# hist(dNdS.Ninteract)
# hist(N_interact$dN)
# hist(N_interact$dS)
# hist(dNdS.Normal)
# hist(Normal$dN)
# hist(Normal$dS)
# dev.off()
# png("F:/Work/Now_work/2021_05_18.ben_/10.dNdS/05.values/Chinese_stats.jpg", width = 1200, height = 800)
# par(mfrow=c(1,3))
# hist(dNdS)
# hist(N_interact$dN)
# hist(N_interact$dS)
# hist(dNdS.Normal)
# hist(Normal$dN)
# hist(Normal$dS)
# dev.off()

data2 <- data[data$dNdS<=5,]

diff.means <- mean(data2$dNdS[data2$type == 1]) - mean(data2$dNdS[data2$type == 0])

one.test <- function(grouping, variable) {
  resampled.group <- sample(grouping)
    mean(variable[resampled.group == 1]) - 
    mean(variable[resampled.group == 0])
}
perm.means <- replicate(1000,one.test(data2$type, data2$dNdS))
sig <- sum(perm.means > diff.means)

