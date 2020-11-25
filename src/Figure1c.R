# ##############################################################################
#
##  Figure 1c Phylum composition
#
# ##############################################################################

L2 <- read.csv('Figure1c.csv',header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE)
###package
#install.packages("reshape2")
#install.packages('plyr')
library(ggplot2)
library(reshape2)
library(plyr)

a <- melt(L2_bar, id.vars = 'Phylum')
L2_bar <- cbind(L2[1],L2[2:12]*100)
L2_bar$Phylum <- reorder(L2_bar$Phylum, L2$`N-FR`)
#write.csv(L2_bar,'L2_bar_plot.csv')
col.hm <- c('#9D9B9C','#550A46','#902044','#CB414B','#D16F7C','#E56B46','#F4A862','#F6DB86','#DFE899',
            '#A4D5A0','#62BA9F','#3681AD','#5A4C97')

g <- ggplot(melt(L2_bar, id.vars = 'Phylum'), aes(x=variable, y=value, fill=Phylum))  +
  geom_bar(stat = "identity", width=0.8, col='black')  +
  #  theme_pander() +
  scale_fill_manual(values = col.hm) +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(hjust=1, angle=90),
        panel.background = element_rect(fill=NULL, colour = 'white')) 

facet_grid(. ~ drv)
ggsave(g, filename = 'figure/Figure1c.pdf', 
       width = 5.5, height = 6, useDingbats=FALSE)
