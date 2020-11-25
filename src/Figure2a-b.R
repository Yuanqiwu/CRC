# ##############################################################################
#
##  Figure 2a-b Biomarker fold change
#
# ##############################################################################

library(ggthemes)
library(ggplot2)
markers <- read.csv('data/figure2a.csv', head = TRUE, 
                    check.names = FALSE, stringsAsFactors = FALSE)
markers_plot <- data.frame(species = markers$species, fold_change = markers$fold_change)
markers_plot$species <- factor(markers_plot$species)
g1 <- ggplot(markers_plot,aes(reorder(species, fold_change), fold_change))+
  geom_bar(aes(fill=factor((fold_change>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#D55E00"), guide=FALSE) +
  labs(x="", y="fold change" ) +
  theme_pander()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust=0),
        panel.background = element_rect(fill=NULL, colour = 'white')
  )
ggsave(g1, filename = 'figure/Figure2a_species.pdf', 
       width = 5.5, height = 6, useDingbats=FALSE)

markers <- read.csv('data/figure2a.csv', head = TRUE, 
                    check.names = FALSE, stringsAsFactors = FALSE)
markers_plot <- data.frame(geneus = markers$geneus, fold_change = markers$fold_change)
markers_plot$geneus <- factor(markers_plot$geneus)
g2 <- ggplot(markers_plot,aes(reorder(geneus, fold_change), fold_change)) +
  geom_bar(aes(fill=factor((fold_change>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#D55E00"), guide=FALSE) +
  labs(x="", y="fold change" ) +
  theme_pander()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust=0),
        panel.background = element_rect(fill=NULL, colour = 'white')
  )
ggsave(g2, filename = 'figure/Figure2a_geneus.pdf', 
       width = 5.5, height = 6, useDingbats=FALSE)
