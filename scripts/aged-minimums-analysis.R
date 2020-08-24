#looking at minimums at different time periods
library(tidyverse)
library(ggpubr)
library(ggthemes)
theme_set(theme_minimal())

mins10k = read.csv("output/results/minimums10k.csv")
mins50k = read.csv("output/results/minimums50k.csv")
mins100k = read.csv("output/results/minimums100k.csv")


reldif10k = mins10k %>% filter(prop != "all")
summary(reldif10k$rel.dif)
reldif50k = mins50k %>% filter(prop != "all")
summary(reldif50k$rel.dif)
reldif100k = mins100k %>% filter(prop != "all")
summary(reldif100k$rel.dif)

reldif = rbind(reldif10k, reldif50k, reldif100k)
reldif$age = ifelse(reldif$prop == "491", "10 ka", 
                    ifelse(reldif$prop == "96", "50 ka", "100 ka"))
reldif$age = factor(reldif$age, levels = c("10 ka", "50 ka", "100 ka"))

reldif2 = rbind(mins10k, mins50k, mins100k)
reldif2$age = ifelse(reldif2$prop == "491", "10 ka", 
                    ifelse(reldif2$prop == "96", "50 ka", 
                           ifelse(reldif2$prop == 12,"100 ka", "all")))
reldif2$age = factor(reldif2$age, levels = c("10 ka", "50 ka", "100 ka", "all"))

with(reldif, pairwise.wilcox.test(rel.dif, age))
with(reldif2, pairwise.t.test(dist.km, age))

testpairs = list(c("50 ka", "100 ka"),
                 c("10 ka", "50 ka"),
                 c("10 ka", "100 ka"))
bplot = ggplot(reldif, aes(x = age, y = dist.km, group = age)) +
  geom_boxplot() +
  stat_compare_means(comparisons = testpairs, method = "wilcox.test", 
                     label = "p.signif") +
  labs(x = "Modeled age of deposits", y = "Distance (km) to first minimum")

ggsave("output/figures/modeled-ages-distdif.png", bplot, 
       dpi = 600, width = 5, height = 4.5)



testpairs = list(c("10 ka", "50 ka"), 
                 c("50 ka", "100 ka"), 
                 c("10 ka", "100 ka"), 
                 c("10 ka", "all"), 
                 c("50 ka", "all"), 
                 c("100 ka", "all"))
bplot2 = ggplot(reldif2, aes(x = age, y = dist.km, group = age)) +
  geom_boxplot() +
  stat_compare_means(comparisons = testpairs, method = "t.test", 
                     label = "p.signif") +
  labs(x = "Modeled age of deposits", y = "Distance to first minimum (km)")
  
cplot = cowplot::plot_grid(bplot2, bplot, labels = c("a)", "b)"))

ggsave("output/figures/modeled-ages-comparison.png", cplot, 
       dpi = 600, width = 11, height = 4.5)
