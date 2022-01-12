####### Summarizing simulation results with tables and graphs

library(tidyverse)
library(ggplot2)
library(here)
library(xtable)

### selection function 1, extreme selection, theta = 0.5
sim1.10 <- readRDS(here("R", "Results", "func1.extreme.big.10.rds"))
sim1.20 <- readRDS(here("R", "Results", "func1.extreme.big.20.rds"))
sim1.40 <- readRDS(here("R", "Results", "func1.extreme.big.40.rds"))
sim1.80 <- readRDS(here("R", "Results", "func1.extreme.big.80.rds"))

### selection function 1, moderate selection, theta = 0.5
sim2.10 <- readRDS(here("R", "Results", "func1.moderate.big.10.rds"))
sim2.20 <- readRDS(here("R", "Results", "func1.moderate.big.20.rds"))
sim2.40 <- readRDS(here("R", "Results", "func1.moderate.big.40.rds"))
sim2.80 <- readRDS(here("R", "Results", "func1.moderate.big.80.rds"))

### selection function 1, extreme selection, theta = 0.1
sim3.10 <- readRDS(here("R", "Results", "func1.extreme.small.10.rds"))
sim3.20 <- readRDS(here("R", "Results", "func1.extreme.small.20.rds"))
sim3.40 <- readRDS(here("R", "Results", "func1.extreme.small.40.rds"))
sim3.80 <- readRDS(here("R", "Results", "func1.extreme.small.80.rds"))

### selection function 1, moderate selection, theta = 0.1
sim4.10 <- readRDS(here("R", "Results", "func1.moderate.small.10.rds"))
sim4.20 <- readRDS(here("R", "Results", "func1.moderate.small.20.rds"))
sim4.40 <- readRDS(here("R", "Results", "func1.moderate.small.40.rds"))
sim4.80 <- readRDS(here("R", "Results", "func1.moderate.small.80.rds"))

### put stacked info into same format as other model info, including iteration number
sim1.10$stacked <- sim1.10$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 10)
sim1.20$stacked <- sim1.20$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 20)
sim1.40$stacked <- sim1.40$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 40)
sim1.80$stacked <- sim1.80$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 80)

sim2.10$stacked <- sim2.10$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 10)
sim2.20$stacked <- sim2.20$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 20)
sim2.40$stacked <- sim2.40$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 40)
sim2.80$stacked <- sim2.80$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 80)

sim3.10$stacked <- sim3.10$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 10)
sim3.20$stacked <- sim3.20$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 20)
sim3.40$stacked <- sim3.40$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 40)
sim3.80$stacked <- sim3.80$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 80)

sim4.10$stacked <- sim4.10$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 10)
sim4.20$stacked <- sim4.20$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 20)
sim4.40$stacked <- sim4.40$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 40)
sim4.80$stacked <- sim4.80$stacked %>%
  mutate(model = "stacked",
         iteration = 1:200,
         theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 80)



# convert columns to correct data type
sim1.10$models <- sim1.10$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 10)
sim1.20$models <- sim1.20$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 20)
sim1.40$models <- sim1.40$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 40)
sim1.80$models <- sim1.80$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Extreme",
         Avg.sample = 80)

sim2.10$models <- sim2.10$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 10)
sim2.20$models <- sim2.20$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 20)
sim2.40$models <- sim2.40$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 40)
sim2.80$models <- sim2.80$models %>%
  type_convert() %>%
  mutate(theta = 0.5,
         Selection = "Moderate",
         Avg.sample = 80)

sim3.10$models <- sim3.10$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 10)
sim3.20$models <- sim3.20$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 20)
sim3.40$models <- sim3.40$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 40)
sim3.80$models <- sim3.80$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Extreme",
         Avg.sample = 80)

sim4.10$models <- sim4.10$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 10)
sim4.20$models <- sim4.20$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 20)
sim4.40$models <- sim4.40$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 40)
sim4.80$models <- sim4.80$models %>%
  type_convert() %>%
  mutate(theta = 0.1,
         Selection = "Moderate",
         Avg.sample = 80)


##### Put all the simulations into a single dataframe
sim.all <- bind_rows(sim1.10$stacked, sim1.10$models, sim1.20$stacked, sim1.20$models,
                     sim1.40$stacked, sim1.40$models, sim1.80$stacked, sim1.80$models,
                     sim2.10$stacked, sim2.10$models, sim2.20$stacked, sim2.20$models,
                     sim2.40$stacked, sim2.40$models, sim2.80$stacked, sim2.80$models,
                     sim3.10$stacked, sim3.10$models, sim3.20$stacked, sim3.20$models,
                     sim3.40$stacked, sim3.40$models, sim3.80$stacked, sim3.80$models,
                     sim4.10$stacked, sim4.10$models, sim4.20$stacked, sim4.20$models,
                     sim4.40$stacked, sim4.40$models, sim4.80$stacked, sim4.80$models) %>%
  arrange(theta, Selection, Avg.sample, iteration)

### create summary stuff for each model/simulation combo
summary.all <- sim.all %>%
                   group_by(theta, Selection, Avg.sample, model) %>%
                   summarize(bias = mean(as.numeric(est.mean)) - theta,
                             avg.sd = mean(as.numeric(sd)),
                             cover = sum(ci.lower < theta & ci.upper > theta) / n(),
                             length = mean(as.numeric(ci.upper) - as.numeric(ci.lower)),
                             rmse = sqrt(mean((as.numeric(est.mean) - theta)^2))) %>%
  arrange(theta, Selection, Avg.sample)

plot.bias <- summary.all %>%
  ggplot(mapping = aes(x = as.factor(Avg.sample), y = abs(bias), group = model, color = model,
                       linetype = (model == "stacked"))) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_wrap( ~ Selection + theta, scales = "free", labeller = label_both) +
  scale_linetype_manual(values = c(1, 3)) +
  scale_color_discrete(labels = c("Bai", "Mavridis", "One-sided (1)", "One-sided (2)",
                                  "One-sided (3)", "One-sided (4)", "Stacked", "Standard",
                                  "Two-sided (1)", "Two-sided (2)"),
                       name = "Model") +
  theme_bw() +
  guides(linetype = "none") +
  labs(x = "Average sample size", y = "|bias|") +
  theme(strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))
  
plot.rmse <- summary.all %>%
  ggplot(mapping = aes(x = as.factor(Avg.sample), y = rmse, group = model, color = model,
                       linetype = (model == "stacked"))) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ Selection + theta, scales = "free", labeller = label_both) +
  scale_linetype_manual(values = c(1, 3)) +
  scale_color_discrete(labels = c("Bai", "Mavridis", "One-sided (1)", "One-sided (2)",
                                  "One-sided (3)", "One-sided (4)", "Stacked", "Standard",
                                  "Two-sided (1)", "Two-sided (2)"),
                       name = "Model") +
  theme_bw() +
  guides(linetype = "none") +
  labs(x = "Average sample size", y = "RMSE") +
  theme(strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

plot.coverage <- summary.all %>%
  ggplot(mapping = aes(x = as.factor(Avg.sample), y = cover, group = model, color = model,
                       linetype = (model == "stacked"))) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ Selection + theta, scales = "free", labeller = label_both) +
  scale_linetype_manual(values = c(1, 3)) +
  scale_color_discrete(labels = c("Bai", "Mavridis", "One-sided (1)", "One-sided (2)",
                                  "One-sided (3)", "One-sided (4)", "Stacked", "Standard",
                                  "Two-sided (1)", "Two-sided (2)"),
                       name = "Model") +
  theme_bw() +
  guides(linetype = "none") +
  labs(x = "Average sample size", y = "95% interval coverage") +
  geom_hline(yintercept = 0.95, linetype = 2, color = "black") +
  theme(strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))


# 
# ggsave(here("Manuscript", "plot.rmse.pdf"), plot = plot.rmse, 
#        width = 6, height = 5, units = "in")
# ggsave(here("Manuscript", "plot.bias.pdf"), plot = plot.bias, 
#        width = 6, height = 5, units = "in")
# ggsave(here("Manuscript", "plot.coverage.pdf"), plot = plot.coverage, 
#        width = 6, height = 5, units = "in")

### Generate latex tables
print(xtable(sim.1.summary[,-7], digits = 3), include.rownames = FALSE)
print(xtable(sim.2.summary[,-7], digits = 3), include.rownames = FALSE)
print(xtable(sim.3.summary[,-7], digits = 3), include.rownames = FALSE)






# graphic for probability of publication in simulation 1

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
  


 
