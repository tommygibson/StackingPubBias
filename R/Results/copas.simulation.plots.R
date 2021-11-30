####### summary graphs and stuff from stacking simulation

library(tidyverse)
library(here)

stacking <- readRDS(here("R", "Results", "small.stack.simulation.rds"))

sim.dat <- stacking$sim.dat
sim.dat$correction <- sim.dat$select.mean - sim.dat$stacked.mean
sim.dat$stack.vs.full <- sim.dat$stacked.mean - sim.dat$full.mean
sim.dat$index <- 1:50

stack.diff <- with(stacking$sim.dat,
                   cbind.data.frame(c(stacked.mean, select.mean, full.mean),
                                    rep(c("stacked", "selected", "full"), each = length(stacked.mean)),
                                    rep(1:length(stacked.mean), 3),
                                    rep(rep(c("1 to 25", "26 to 50"), each = 25), 3)))
names(stack.diff) <- c("Mean", "source", "index", "set")

stack.diff %>%
  ggplot() +
  geom_boxplot(aes(x = Mean, fill = source)) +
  facet_wrap( ~ source, ncol = 1) +
  labs(title = "Distribution of data vs stacked estimates")

stack.diff %>%
  ggplot(aes(x = Mean, y = index, color = source)) +
  geom_point() 

sim.dat %>%
  ggplot(aes(x = correction)) +
  geom_histogram()

sim.dat %>%
  ggplot(aes(x = stack.vs.full, y = index)) +
  geom_point()
  
