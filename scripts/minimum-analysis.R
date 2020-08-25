library(tidyverse)
library(bbmle)
library(mgcv)
library(stargazer)
library(rcompanion)
theme_set(theme_minimal())

####Total Dataset Analysis####
mins1 = read.csv("output/results/minimums.csv")
mins2 = read.csv("output/results/small_minimums.csv")

### regression on datasets separately
fit100 = lm(log(dist.km) ~ prop, data = mins1)
plot(fit100, which = 2)
summary(fit100)
#mins2_2 = mins2 %>% filter(prop != 100)
fit10 = lm(log(dist.km) ~ prop, data = mins2)
plot(fit10, which = 2)
summary(fit10)
stargazer(fit100, fit10, type = "text")

pd100 = data.frame(prop = seq(0, 100))
pd100$dist.km = exp(predict(fit100, newdata = pd100))
#pd100$dist.km = predict(fit100, newdata = pd100)
# plot100 = ggplot(data = mins1, aes(x = prop, y = log(dist.km))) +
#   geom_boxplot(aes(group = prop), alpha = 0.5) +
#   geom_line(data = pd100, aes(x = prop, y = dist.km), color = "green", size = 1) +
#   labs(x = "Percentage of total points", y = "log(Distance to first minimum (km))")
# plot(plot100)

pd10 = data.frame(prop = seq(0, 100))
pd10$dist.km = exp(predict(fit10, newdata = pd10))
#pd10$dist.km = predict(fit10, newdata = pd10)
# plot10 = ggplot(data = mins2, aes(x = prop, y = log(dist.km))) +
#   geom_boxplot(aes(group = prop), alpha = 0.5) +
#   geom_line(data = pd10, aes(x = prop, y = dist.km), color = "green", size = 1) +
#   labs(x = "Percentage of total points", y = "log(Distance to first minimum (km))")
# plot(plot10)


mins = rbind(mins1, mins2)
plotboth = ggplot() +
  geom_boxplot(data = mins, aes(x = prop, y = dist.km, group = prop), alpha = 0.25) +
  geom_line(data = pd10, aes(x = prop, y = dist.km), color = "blue", size = 0.5) +
  geom_line(data = pd100, aes(x = prop, y = dist.km), color = "red", size = 0.5) +
  labs(x = "Percentage of total points", y = "Distance to first minimum (km)")
plot(plotboth)

# with(mins, plot(dist.km ~ prop))
# with(mins, plot(sim ~ prop))
# with(mins, plot(rel.dif ~ prop))

expfit1 = nls(dist.km ~ a*(prop^b), start = c(a = 1, b = 1), data = mins)
summary(expfit1)
newmins$p2 = 1322*(newmins$prop^-0.04148)
summary(lm(dist.km ~ p2, data = newmins))

epdf = data.frame(prop = seq(1, 100))
epdf$dist.km = predict(expfit1, newdata = epdf)
eplot1 = ggplot(data = mins, aes(x = prop, y = dist.km)) +
  geom_boxplot(aes(group = prop), alpha = 0.5, size = 0.5) +
  geom_line(data = epdf, aes(x = prop, y = dist.km), color = "green", size = 0.5) +
  labs(x = "Percentage of total points", y = "Distance to first minimum (km)")
plot(eplot1)

cowplot::plot_grid(eplot1, plotboth, labels = c("a)", "b)"), ncol = 1)
# 
# ggsave("output/figures/absolute-dist_regs.png", cowplot::plot_grid(eplot1, plotboth, labels = c("a)", "b)"), ncol = 1), 
#        dpi = 600, width = 5, height = 7)

# expfit2 = nls(rel.dif ~ a*(prop^b), start = c(a = 1, b = 1), data = newmins)
# summary(expfit2)
# newmins$p3 = 170.10374*(newmins$prop^-0.53321)
# summary(lm(rel.dif ~ p3, data = newmins))
# 
# epdf2 = data.frame(prop = seq(min(mins$prop), max(mins$prop)))
# epdf2$rel.dif = predict(expfit2, newdata = epdf2)
# eplot2 = ggplot(data = newmins, aes(x = prop, y = rel.dif)) +
#   geom_boxplot(aes(group = prop), alpha = 0.5, size = 0.5) +
#   geom_line(data = epdf2, aes(x = prop, y = rel.dif), color = "green", size = 0.5) +
#   labs(x = "Percentage of total points", y = "Relative difference in distance (km) to first minimum")

# fit1 = lm(dist.km ~ prop, data = mins)
# plot(fit1, which=2)
# summary(fit1)
# 
# fit2 = lm(sim ~ prop, data = mins)
# plot(fit2, which = 2)
# plot(fit2, which = 1)
# summary(fit2)
# 
# fit3 = lm(rel.dif ~ prop, data = mins)
# plot(fit3, which = 2)
# plot(fit3, which = 1)
# summary(fit3)

####Germany Analysis####
de_mins1 = read.csv("output/results/minimums_de.csv")
de_mins2 = read.csv("output/results/small_minimums_de.csv")
de_mins = rbind(de_mins1, de_mins2)

fit4 = lm(dist.km ~ prop, data = de_mins)
plot(fit4, which = 2)
summary(fit4)
fit5 = lm(log(dist.km) ~ prop, data = de_mins)
plot(fit5, which = 2)
summary(fit5)
stargazer(fit5, title = "Linear Regression Results", type = "text")

ppdf = data.frame(prop = seq(0, 100))
ppdf$dist.km = exp(predict(fit5, newdata = ppdf, type = "response"))
#ppdf$dist.km = predict(fit5, newdata = ppdf, type = "response")
plot1 = ggplot(data = de_mins, aes(x = prop, y = dist.km)) +
  geom_boxplot(aes(group = prop), alpha = 0.5) +
  geom_line(data = ppdf, aes(x = prop, y = dist.km), color = "green", size = 0.75) +
  labs(x = "Percentage of total points", y = "Distance to first minimum (km)")
plot(plot1)
# ggsave("output/figures/distance-to-min_de.png", plot1, 
#        dpi = 600, width = 5, height = 4)


####unused analysis####
# ppdf2 = data.frame(prop = ppdf$prop)
# ppdf2$rel.dif = predict(fit4, newdata = ppdf2)
# plot2 = ggplot(data = newmins, aes(x = prop, y = rel.dif)) +
#   geom_boxplot(aes(group = prop), alpha = 0.5, size = 1) +
#   geom_line(data = ppdf2, aes(x = prop, y = rel.dif), color = "green", size = 1) +
#   labs(x = "Percentage of total points", y = "Relative difference in distance (km) to first minimum")
# #ggsave("output/figures/rel-diff-to-min_all-perc.png", plot2)
# 
# cplot = cowplot::plot_grid(plot1, plot2, labels = c("a)", "b)"), ncol = 1, nrow = 2, hjust = -0.25, scale = 0.9)
# ggsave("output/figures/germany-results.png", 
#        dpi = 600, width = 4, height = 6)
# 
# large1 = lm(dist.km ~ prop, data = mins1)
# summary(large1)
# large2 = lm(rel.dif ~ prop, data = mins1)
# summary(large2)
# 
# 
# small1 = lm(dist.km ~ prop, data = mins2)
# summary(small1)
# small2 = lm(rel.dif ~ prop, data = mins2)
# summary(small2)
# 
# 
# plotmins = mins %>% filter(prop == 100 | prop == 50 | prop == 10 |
#                              prop == 5 | prop == 1)
# with(plotmins, pairwise.t.test(dist.km, prop))
# ggplot(plotmins, aes(x = as.factor(prop), y = dist.km, group = as.factor(prop))) +
#   geom_boxplot()
# 
# stargazer(fit4, fit5, type = "html",
#           dep.var.labels = c("Relative difference in distance",
#                              "Absolute distance"),
#           covariate.labels = c("Percentage of points"),
#           out = "output/models.htm")
# 
# 
# 
# 
# 
# # ggsave("output/figures/power-fit-rel-dist.png", eplot2, 
# #        width = 7, height = 4, dpi = 600)
# 
# cor(fitted(expfit2), epdf2$rel.dif)
# fitted(expfit2)

