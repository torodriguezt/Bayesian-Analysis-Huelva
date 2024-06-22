listcombo <- unlist(sapply(0:7, function(x) combn(7, x, simplify=FALSE)), recursive=FALSE)
predterms <- lapply(listcombo, function(x) paste(c("b0", "f(s, model = spde)",c("ind","mar","agric","refi", "urban", "phospho", "bare")[x]),collapse="+"))

coefm <- matrix(NA, 128, 3)

for(i in 1:128){
  formula <- as.formula(paste("y ~ ", predterms[[i]]))
  result <- inla(formula, family = "lognormal", 
                 data = inla.stack.data(stk.full),
                 control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
                 control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                        return.marginals.predictor = TRUE,
                                        config = TRUE),
                 control.fixed=list(mean = 0,prec = 0.01,
                                    mean.intercept = 0,prec.intercept = 0.01))
  coefm[i,1] <- result$dic$dic
  coefm[i,2] <- result$waic$waic
  coefm[i,3] <- -sum(log(result$cpo$cpo[1:138]))
}

rownames(coefm) <- predterms
colnames(coefm) <- c("DIC", "WAIC", "CPO")
round(coefm, 3)
