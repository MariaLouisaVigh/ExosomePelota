rm(list=ls(all=TRUE))

library(openxlsx)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggiraphExtra)
library(ggiraph)
library(lme4)
library(report)
library(svglite)
library(gof)
library(emmeans)


load(file = "Bindingenergy.RDat")

#make a subset having only one binding energy per target and summing the miRNA RPM of these

df6 <-  df3%>%
  group_by(target_id)%>%
  mutate(count = n_distinct(mseq5)>1)%>% 
  filter(., count=="FALSE")%>%
  mutate(lsumRPM5 = log(sum(miRNA_RPM)))%>%
  as.data.frame(.,)  

# Make dataset with only miRNA families that has targets making siRNA and not
dif <-df3 %>% 
  dplyr::group_by(miRNAfamily, sign.diff.target) %>%
  dplyr::summarise(mean_Gibbs5a = mean(Gibbs5a, na.rm=TRUE))%>%
  tidyr::spread(., key= "sign.diff.target", value = "mean_Gibbs5a" )%>%
  as.data.frame(.)

dif2 <- dif[(!is.na(dif$no)&(!is.na(dif$yes))),]
dif3 <- merge(df3, dif2, by="miRNAfamily")

#Plot ordered miRNA expression and targets making siRNA or not
A <- ggplot(data =df6, aes(x = reorder(mt,lsumRPM5) , y = lsumRPM5, fill = sign.diff.target)) +
  geom_bar(stat = "identity") +
  ggtitle("miRNA and targets ordered by miRNA logRPM")+
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size=10),
        axis.text.y = element_text(size=10),legend.title=element_text(size=10),
        legend.text=element_text(size=10), axis.text.x = element_blank(), 
        axis.title.y = element_text(size=10, vjust = 1), plot.margin=unit(c(1,1,1,1),"cm"))
A


B <- ggplot(data =df6, aes(x = reorder(mt,(lsumRPM5*-Gibbs5a)) , y = (lsumRPM5*-Gibbs5a), fill = sign.diff.target)) +
  geom_bar(stat = "identity") +
  ggtitle("miRNA and targets ordered by miRNA logRPM * -deltaG")+
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(size = 10),axis.title=element_text(size=10),
        axis.text.y = element_text(size=10),legend.title=element_text(size=10),
        legend.text=element_text(size=10), axis.text.x = element_blank(), 
  axis.title.y = element_text(size=10, vjust = 1), plot.margin=unit(c(1,1,1,1),"cm"))
B

# Logistic models of the effect of miRNA expression on siRNA production 

m1 <- glm(siRNA ~ lsumRPM5 , data = df6, family =binomial)

plot(cumres(m1))
summary(m1)
drop1(m1, test = "Chisq")
report(m1)

# Saturated logistic models of the effect of miRNA expression and 
#Gibbs5a on siRNA production 
m2 <- glm(siRNA ~ lsumRPM5 * Gibbs5a , data = df6, family =binomial)

plot(cumres(m2))
summary(m2)
drop1(m2, test = "Chisq")
report(m2)

# Logistic models of the effect of miRNA expression and Gibbs5a
# as covariate on siRNA production 
m3 <- glm(siRNA ~ I(lsumRPM5* Gibbs5a) , data = df6, family =binomial)

plot(cumres(m3))
summary(m3)
drop1(m3, test = "Chisq")
report(m3)

# Linear mixed model to test the effect of siRNA or not on 
# binding energy within miRNA family
m4 <- lmer(Gibbs5a ~sign.diff.target + (1|miRNAfamily), data = dif3)

qqnorm(resid(m4))
plot(m4)
summary(m4)
drop1(m4, test = "Chisq")
report(m4)

#Plot binding energy within miRNA families
C <- ggplot(data=dif3, aes(sign.diff.target, -Gibbs5a, color =sign.diff.target))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  scale_y_continuous(position = "right") +
  geom_boxplot(fill="transparent", size=0.8, outlier.shape=1, outlier.size=1, width=0.5, position = position_dodge(), coef=10) +
  geom_point(shape=21, col="black", size = 2,  fill="transparent", position=position_dodge2(width=0.5), alpha=.5) +
  theme(strip.text = element_text(size = 10))+
  cowplot::theme_cowplot() +
  facet_wrap(~miRNAfamily)

C
