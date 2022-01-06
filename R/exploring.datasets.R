### picking some datasets for illustrating publication bias

library(metafor)
library(tidyverse)

plot(dat.begg1989$sei ~ dat.begg1989$yi)
plot(sqrt(dat.knapp2017$vi) ~ dat.knapp2017$yi)
plot(dat.lopez2019$se ~ dat.lopez2019$diff)
dat.anand1999 %>%
  mutate(si = sqrt(vi)) %>%
  ggplot(mapping = aes(x = yi, y = si, color = as.factor(pubstatus))) +
  geom_point()

plot(sqrt(dat.bangertdrowns2004$vi) ~ dat.bangertdrowns2004$yi) # better for meta-regression

dat.landenberger2005 %>%
  mutate(log.or = log((n.cbt.non * n.ctrl.rec) / (n.cbt.rec * n.ctrl.non)),
         se.log.or = sqrt(1 / n.cbt.non + 1 / n.ctrl.rec + 1 / n.cbt.rec + 1 / n.ctrl.non)) %>%
  ggplot(mapping = aes(x = log.or, y = se.log.or)) + 
  geom_point()
  
  
  
dat.3 <- dat.bornmann2007 %>%
  mutate(a = maward,
         b = waward,
         c = mtotal - maward,
         d = wtotal - waward,
         nm = maward,
         nw = waward,
         log.or = log((a * d) / (b * c)),
         se.log.or = sqrt(1 / a + 1 / b + 1 / c + 1 / d))
plot(dat.born$se.log.or ~ dat.born$log.or)
