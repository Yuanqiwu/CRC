# ##############################################################################
#
##  Figure 1b beta diversity
#
# ##############################################################################
# Packages
library("labdsv")
library("coin")
library("vegan")
library("yaml")
library("ggpubr")
library("cowplot")
library("tidyverse")
#import data
df.plot <- read.csv(file = 'data/Figure1b.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                                check.name = FALSE)
#main plot
g.main <- df.plot %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape= Group ,col= Study)) +
  geom_point(size = 1.5 )  +
  xlab(axis.1.title) + ylab(axis.2.title) +
  scale_colour_manual(values=c('CA-CRC' = '#085F63', 'FR-CRC' = '#49BEB7',
                               'US1-CRC' = '#FF5959', 'US2-CRC' = '#FACF5A'),guide = FALSE) +
  scale_shape_manual(values=c(15, 16, 17),guide=FALSE) +
  scale_x_continuous(position='top') +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks=element_blank(), axis.text = element_blank(),axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank())
#legend.title=element_text(size =15),legend.text=element_text(size = 15),

# study boxplot axis 1
g.s.1 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(c('CA-CRC' = '#085F63', 'FR-CRC' = '#49BEB7',
                                            'US1-CRC' = '#FF5959', 'US2-CRC' = '#FACF5A')))) %>% 
  ggplot(aes(y=Axis1, x=Study, fill=Study)) +
  xlab(paste0('Study\n','')) +
  geom_boxplot() +
  scale_fill_manual(values=c('CA-CRC' = '#085F63', 'FR-CRC' = '#49BEB7',
                             'US1-CRC' = '#FF5959', 'US2-CRC' = '#FACF5A'), guide=FALSE) +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y  = element_text(size = 10),
        panel.grid = element_blank()) + 
  coord_flip()

#stat_compare_means(label = "p.format", method = "kruskal.test", size=2) +
# study boxplot axis 2
g.s.2 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(c('CA-CRC' = '#085F63', 'FR-CRC' = '#49BEB7',
                                            'US1-CRC' = '#FF5959', 'US2-CRC' = '#FACF5A')))) %>% 
  ggplot(aes(y=Axis2, x=Study, fill=Study)) + 
  xlab(paste0('Group\n','')) +
  geom_boxplot() + 
  scale_fill_manual(values=c('CA-CRC' = '#085F63', 'FR-CRC' = '#49BEB7',
                             'US1-CRC' = '#FF5959', 'US2-CRC' = '#FACF5A'), guide = FALSE) +
  
  scale_x_discrete(position='top') +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x  = element_text(size = 10),
        panel.grid = element_blank())

# group plot axis1
g.g.1 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis1, fill=Group)) +
  xlab(paste0('Study\n','')) +
  geom_boxplot() +
  scale_fill_manual(values= c('Normal' = '#C2C2C2', 'Adenoma' = '#226597', 'Cancer' = '#d9544d'),guide= FALSE) + 
  ylab(axis.1.title) +
  
  theme(axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y=element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x  = element_blank(),
        legend.title=element_text(size =15),legend.text=element_text(size = 15),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()

# group plot axis2
g.g.2 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis2, fill=Group)) +
  xlab(paste0('Group\n','')) +
  geom_boxplot() +
  scale_fill_manual(values=c('Normal' = '#C2C2C2', 'Adenoma' = '#226597', 'Cancer' = '#d9544d'), guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  
  ylab(axis.2.title) + 
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())

pdf('Figures/Figure1b.pdf', useDingbats = FALSE)
a <- plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL,NULL,g.g.1, NULL, NULL,
               nrow=3,
               rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))
dev.off()