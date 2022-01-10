##### Functions for StackingPubBias

summary.func <- function(x){
  c(mean(x), sd(x), quantile(x, c(.025, .975)))
}

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
