#' @export
RDA_factors <- function(names){
  labs = as.data.frame(do.call("rbind", strsplit(names,"_")))
  colnames(labs) = lapply(strsplit(comparisons,"\\\\"),function(x){
    paste(x,collapse="-")
  })
  mod0 <- rda(t(expr.matrix) ~ 1, labs)
  mod1 <- rda(t(expr.matrix) ~ ., labs)
  #mod <- ordistep(mod0, scope=formula(mod1))
  #plot(mod)
  step.res <- ordiR2step(mod0, mod1, perm.max = 200)
  step.res$anova  # Summary table
  plot(step.res)
}
