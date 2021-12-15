####### Summarizing simulation results with tables and graphs

library(tidyverse)
library(ggplot2)
library(here)
library(xtable)

sim1 <- readRDS(here("R", "Results", "sim.results.1.rds"))
sim2 <- readRDS(here("R", "Results", "sim.results.2.rds"))

sim1$stacked <- sim1$stacked %>%
  mutate(model = "stacked",
         iteration = 1:1000,
         simulation = 1)
sim2$stacked <- sim2$stacked %>%
  mutate(model = "stacked",
         iteration = 1:1000,
         simulation = 2)

sim1$models <- sim1$models %>%
  mutate(simulation = 1) %>%
  type_convert()
  mutate(simulation = 1)
sim2$models <- sim2$models %>%
  mutate(simulation = 2) %>%
  type_convert() 



all.models.1 <- bind_rows(sim1$stacked, sim1$models, sim2$stacked, sim2$models) %>%
  arrange(simulation, iteration)

baddies.1 <- unique(all.models.1$iteration[is.na(all.models.1$est.mean) & all.models.1$simulation==1])
baddies.2 <- unique(all.models.1$iteration[is.na(all.models.1$est.mean) & all.models.1$simulation==2])

sim.1.summary <- all.models.1 %>%
                   filter(!(iteration %in% baddies.1),
                          simulation == 1) %>%
                   group_by(model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - sim1$theta0,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < 0.4 & ci.upper > 0.4) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)),
                             base.studies = 30)
rownames(sim.1.summary) <- NULL
sim.2.summary <- all.models.1 %>%
                   filter(!(iteration %in% baddies.2),
                          simulation == 2) %>%
                   group_by(model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - sim1$theta0,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < 0.4 & ci.upper > 0.4) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)),
                             base.studies = 60)

sim.full.summary <- bind_rows(sim.1.summary, sim.2.summary)

sim.full.summary %>%
  ggplot(aes(x = base.studies, y = rmse, color = model)) +
  geom_point() + 
  geom_line()

sim.full.summary %>%
  ggplot(aes(x = base.studies, y = avg.sd, color = model)) +
  geom_point() + 
  geom_line()

### Generate latex tables
# print(xtable(sim.1.summary[,-7], digits = 3), include.rownames = FALSE)
# print(xtable(sim.2.summary[,-7], digits = 3), include.rownames = FALSE)

 
