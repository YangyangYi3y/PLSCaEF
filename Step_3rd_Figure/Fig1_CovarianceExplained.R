
library(R.matlab)
library(ggplot2)

Folder = 'D:/projects/ABCD/EF/PLS/results/ef_task23/combat/2fold';

# All items
WeightMat = readMat(paste0(Folder, '/Comp_corr.mat'));

CovarianceExplained = WeightMat$CovarianceExplained.Median;
CovarianceExplained = CovarianceExplained / sum(CovarianceExplained);
data = data.frame(ComponentNumber = as.numeric(c(1:23)));
data$CovExplained_All = as.numeric(CovarianceExplained);
component1 <- data[data$ComponentNumber=='1',]
component2 <- data[data$ComponentNumber=='2',]

dev.off()
dev.new()
windowsFonts(arial = windowsFont(family = "Arial"))          

par(pin = c(6,5))
Fig <- ggplot(data, aes(x=ComponentNumber, y=CovExplained_All)) +
    geom_point(shape = 21, stroke = 2, col='grey',fill = 'grey',size = 4) + 
    theme_classic() + 
    labs(x = "Component", y = "Covariance Explained") + 
    theme(axis.title.x = element_text(margin=margin(t=10),face="bold", color='black', size=18,family = "arial"),
        axis.text.x = element_text(margin=margin(t=0), size=8, angle=0))+ 
    theme(axis.title.y = element_text(margin=margin(r=20),angle=90, face="bold", size=18,color='black',family = "arial"),
        axis.text.y = element_text(margin=margin(r =0),size=16)) +
    scale_x_continuous(breaks = seq(1, 23, by=1))+
    scale_y_continuous(breaks = seq(0, 0.4, by=0.02))+
    geom_point(data=component1, 
               shape = 21, stroke = 2,col='#4472c4',fill = 'grey',size = 4)+
    geom_point(data=component2, 
               shape = 21, stroke = 2, col='#ed7d31',fill = 'grey',size = 4)

Fig

ggsave('/CovarianceExplained.tiff', width = 12, height = 10, dpi = 600, units = "cm")

