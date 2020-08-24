library(tidyverse)
theme_set(theme_minimal())

c14 = readxl::read_xlsx("data/radiocarbon-palaeolithic-europe-database-v25-extract.xlsx", skip = 7) %>%
  filter(reliability == "reliable" | is.na(reliability))
dates = data.frame(age = c14$age) %>% filter(!is.na(age)) %>% 
  filter(age <= 900000 & age >= 7000)
summary(dates$age)

min = min(dates$age)
max = max(dates$age)
bins = seq(min, max, by = 500)
dates$bins = cut(dates$age, breaks = bins, labels = bins[-1], include.lowest = T)
dates$binsN = as.numeric(as.character(dates$bins))
freq = dates %>% group_by(binsN) %>% summarize(n = n())

theta = min(freq$n) * 0.5
fit1 = lm(log(n - theta) ~ binsN, data = freq)
alpha = exp(coef(fit1)[1])
beta = coef(fit1)[2]
start = list(alpha = alpha, beta = beta, theta = theta)
model = nls(n ~ alpha * exp(beta * binsN) + theta, data = freq, start = start)
summary(model)

freq$predictor = 412.5*exp(freq$binsN*(-4.073*10^-5))
summary(lm(n ~ predictor, data = freq))

ppdf = data.frame(binsN = seq(min(freq$binsN), max(freq$binsN)))
ppdf$n = predict(model, newdata = ppdf)

ggplot(data = freq, aes(x = binsN, y = n)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_line(data = ppdf, aes(binsN, n), color = "blue")


model2 = nls(n ~ a*((binsN)^c), start = list(a = 0.5, c = 0.5), data = freq)
summary(model2)
freq$p2 = (3.736*10^6)*(freq$binsN^(-1.031))
summary(lm(n ~ p2, data = freq))

ppdf$n2 = predict(model2, newdata = ppdf)
ggplot(data = freq, aes(x = binsN, y = n)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_line(data = ppdf, aes(binsN, n2), color = "blue")


