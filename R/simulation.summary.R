####### Summarizing simulation results with tables and graphs

library(tidyverse)
library(ggplot2)
library(here)
library(xtable)

sim1 <- readRDS(here("R", "Results", "sim.results.1.rds"))
sim2 <- readRDS(here("R", "Results", "sim.results.2.rds"))
sim3 <- readRDS(here("R", "Results", "sim.results.3.rds"))

sim1$stacked <- sim1$stacked %>%
  mutate(model = "stacked",
         iteration = 1:1000,
         simulation = 1)
sim2$stacked <- sim2$stacked %>%
  mutate(model = "stacked",
         iteration = 1:1000,
         simulation = 2)

sim3$stacked <- sim3$stacked %>%
  mutate(model = "stacked",
         iteration = 1:500,
         simulation = 3)

# convert columns to correct data type
sim1$models <- sim1$models %>%
  type_convert()
sim2$models <- sim2$models %>%
  type_convert() 
sim3$models <- sim3$models %>%
  type_convert()



all.models.1 <- bind_rows(sim1$stacked, sim1$models) %>%
  arrange(simulation, iteration)
all.models.2 <- bind_rows(sim2$stacked, sim2$models) %>%
  arrange(simulation, iteration)
all.models.3 <- bind_rows(sim3$stacked, sim3$models) %>%
  arrange(simulation, iteration)

# this is from when we were throwing out iterations where pareto k's were big
# baddies.1 <- unique(all.models.1$iteration[is.na(all.models.1$est.mean) & all.models.1$simulation==1])
# baddies.2 <- unique(all.models.1$iteration[is.na(all.models.1$est.mean) & all.models.1$simulation==2])

sim.1.summary <- all.models.1 %>%
                   group_by(model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - sim1$theta0,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < 0.4 & ci.upper > 0.4) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)),
                             base.studies = 30)

sim.2.summary <- all.models.2 %>%
                   group_by(model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - sim1$theta0,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < 0.4 & ci.upper > 0.4) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)),
                             base.studies = 60)

sim.3.summary <- all.models.3 %>%
                   group_by(model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - sim1$theta0,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < 0.4 & ci.upper > 0.4) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - 0.4)^2)),
                             base.studies = 100)

sim.full.summary <- bind_rows(sim.1.summary, sim.2.summary, sim.3.summary)

sim.full.summary %>%
  ggplot(aes(x = as.factor(base.studies), y = bias, group = model, color = model)) +
  geom_point() + 
  geom_line()

sim.full.summary %>%
  ggplot(aes(x = base.studies, y = avg.sd, color = model)) +
  geom_point() + 
  geom_line()

### Generate latex tables
# print(xtable(sim.1.summary[,-7], digits = 3), include.rownames = FALSE)
# print(xtable(sim.2.summary[,-7], digits = 3), include.rownames = FALSE)


# graphic for probabilty of publication in simulation 1

propensity.func.1 <- function(p){
  p.prop <- ifelse(p < 0.005, 1, 
                   ifelse(p < 0.2, exp(-2 * p), 
                          ifelse(p < 0.5, exp(-4 * p), .1)))
  
  return(p.prop)
}
propensity.func.2 <- function(p){
  p.prop <- ifelse(p < 0.005, 1,
                   ifelse(p < 0.2, exp(-0.5 * p),
                          ifelse(p < 0.5, exp(-1 * p), .5)))
  return(p.prop)
}

y <- seq(-.5, 1, .05)
s <- seq(.2, .8, length.out = length(y))
heatmap.dat <- as.data.frame(expand.grid(y, s)) %>%
                 rename("y" = "Var1",
                        "s" = "Var2") %>%
                 mutate(p = 1 - pnorm(y / s),
                        p.select.1 = as.factor(round(propensity.func.1(p), 2)),
                        p.select.2 = as.factor(round(propensity.func.2(p), 2)))
P1 <- nlevels(heatmap.dat$p.select.1)
P2 <- nlevels(heatmap.dat$p.select.2)

colors1 <- colorRampPalette(c("firebrick", "gold2", "olivedrab3", "royalblue"))(P1)
colors2 <-colorRampPalette(c("gold2", "olivedrab3", "royalblue"))(P2)
(heatmap1 <- heatmap.dat %>%
  ggplot(aes(x = y, y = s, fill = p.select.1)) +
  geom_tile() +
  scale_fill_manual(values = colors1, breaks = levels(heatmap.dat$p.select.1)[seq(1, P1, 7)]) +
  theme_bw())
(heatmap2 <- heatmap.dat %>%
    ggplot(aes(x = y, y = s, fill = p.select.2)) +
    geom_tile() +
    scale_fill_manual(values = colors2, breaks = levels(heatmap.dat$p.select.2)[seq(1, P2, 4)]) +
    theme_bw())

p <- seq(0, 1, 0.005)
p.select.1 <- propensity.func.1(p)
p.select.2 <- propensity.func.2(p)

propensity.dat <- data.frame(cbind(rep(p, 2),
                                   c(p.select.1, p.select.2),
                                   c(rep("Extreme", length(p)), rep("Moderate", length(p)))))
names(propensity.dat) <- c("p.value", "prob.select", "Severity")
(regmap <- propensity.dat %>%
    group_by(Severity) %>%
    type_convert() %>%
    ggplot(aes(x = p.value, y = prob.select, color = Severity)) +
      geom_line() +
      labs(title = "Selection probabilities for Selection Mechanism 1",
           x = "p-value", y = "Selection probability") +
      theme_bw())
ggsave(filename = "SM1.pdf", plot = regmap, path = here("R", "Figures"),
       width = 6, height = 4)
  


 
