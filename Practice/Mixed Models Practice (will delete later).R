
library(lme4)

Estuaries <- read.csv("Estuaries.csv", header = T)

ft.estu <- lmer(Total ~ Modification + (1|Estuary),data = Estuaries, REML=T)

qqnorm(residuals(ft.estu))
