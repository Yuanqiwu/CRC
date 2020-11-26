# ##############################################################################
#
##  Figure 1a confounder
#
# ##############################################################################

#  variance explained by confounders
ss.disease <- apply(CRC_relative, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=metadata %>% pull(Group))

# calculate trimmed mean abundance
t.mean <- apply(CRC_relative, 1, mean, trim=0.1)

# plot study
alpha.meta <- 0.05
CRC_plot.all <- read.csv('data/Figure1a',row.names = 1,header=TRUE,
                            stringsAsFactors = FALSE,check.names = FALSE)
df.plot.study <- CRC_plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  filter(type=='Study')

g2 <- df.plot.study %>%
  ggplot(aes(x=disease, y=meta)) +
  geom_point(aes(size=t.mean, fill=meta.significance), shape=21,
             col=alpha(c('black'), alpha=0.4)) +
  xlab(paste0('Variance explained by Disease\n','species',' average: ',
              formatC(mean(df.plot.study$disease)*100, digits=2), '%')) +
  ylab(paste0('Variance explained by Study\n','species',' average: ',
              formatC(mean(df.plot.study$meta)*100, digits=2), '%')) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(from=0, to=0.1, by=0.05)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_fill_manual(values = alpha(c('grey', '#CC071E'),
                                   alpha=c(0.4, .8)),
                    name=paste0('Significance\n(', alpha.meta,')')) +
  scale_size_area(name='Trimmed mean abundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')

ggsave(g2, filename = 'figure/Figure1a.pdf',
       width = 6, height = 6)

cat('Successfully computed confounder effects in',
    proc.time()[1]-start.time, 'second...\n')

# ##############################################################################
