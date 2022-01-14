### quick tests after simulations
library(tidyverse)

sim4.40$models %>%
  group_by(model) %>%
  summarize(m = mean(as.numeric(est.mean)),
            s = sd(as.numeric(est.mean)),
            ms = mean(as.numeric(sd)),
            cover = sum(ci.lower < 0.1 & ci.upper > 0.1),
            length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
            rmse = sqrt(mean((as.numeric(est.mean) - 0.1)^2)))

