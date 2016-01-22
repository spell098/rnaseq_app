#' RDA plots
#' Creates a biplot
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' expr.matrix <- readRDS('data/expr_matrix_LGVD.rds')
#' RDA_factors(expr.matrix[1:100,])
#' @export
RDA_factors <- function(expr.matrix){
  names <- get_names(expr.matrix)
  labs = as.data.frame(do.call("rbind", strsplit(names,"_")))
  colnames(labs) = apply(labs,2,function(x){
    paste(unique(x),collapse="_")
  })
  mod0 <- rda(t(expr.matrix) ~ 1, labs)
  mod1 <- rda(t(expr.matrix) ~ ., labs)
  mod <- ordistep(mod0, scope=formula(mod1))
  plot(mod)
  step.res <- ordiR2step(mod0, mod1, perm.max = 200)
  step.res$anova  # Summary table
  #plot(step.res)
}
