get.result.vegan <- function(objVegan){
            nbrElements = length(objVegan$colsum)
            Total_variance   =  objVegan[[5]]
            Constraint_Proportion     =  sum( objVegan$CCA$eig)/objVegan[[5]]
            Unconstraint_Proportion   =  sum( objVegan$CA$eig) /objVegan[[5]]
            RDA1                      =  objVegan$CCA$eig[1]   /objVegan[[5]]
            RDA2                      =  objVegan$CCA$eig[2]   /objVegan[[5]]
            RDA3                      =  objVegan$CCA$eig[3]   /objVegan[[5]]
            PC1                       =  objVegan$CA$eig[1]    /objVegan[[5]]
            PC2                       =  objVegan$CA$eig[2]    /objVegan[[5]]
            PC3                       =  objVegan$CA$eig[3]    /objVegan[[5]]
            R2                        <-RsquareAdj(objVegan)$r.squared     # unadjusted Square R
            R2adj                     <-RsquareAdj(objVegan)$adj.r.squared # adjusted Square R
            RDA1_compared_Explained_variance                     =  (RDA1/Constraint_Proportion)*R2adj
            RDA2_compared_Explained_variance                      =  (RDA2/Constraint_Proportion)*R2adj
            RDA3_compared_Explained_variance                      =  (RDA3/Constraint_Proportion)*R2adj
            rs_anova     =   anova.cca(objVegan, step= 1000)
            Pvalue_model =   rs_anova$"Pr(>F)"[1]

            rs_dataframe <- data.frame(nbrElements= nbrElements                            ,
                                      Total_variance   =  Total_variance                   ,
                                      Constrained_Proportion     =  Constraint_Proportion  ,
                                      Unconstrained_Proportion   =  Unconstraint_Proportion,
                                      RDA1                      =  RDA1                     ,
                                      RDA2                      =  RDA2                     ,
                                      RDA3                      =  RDA3                     ,
                                      PC1                       =  PC1                      ,
                                      PC2                       =  PC2                      ,
                                      PC3                       =  PC3                      ,
                                      R2                        = R2                       ,
                                      R2adj                     = R2adj                    ,
                                      RDA1_vs_Explained_variance = RDA1_compared_Explained_variance,
                                      RDA2_vs_Explained_variance = RDA2_compared_Explained_variance,
                                      RDA3_vs_Explained_variance = RDA3_compared_Explained_variance,
                                      Pvalue_model =   Pvalue_model
                                     )
           rs_dataframe = t(rs_dataframe)
           colnames(rs_dataframe) <- "values"
           return(round(rs_dataframe, 3))
}#end of get.vegan.result
