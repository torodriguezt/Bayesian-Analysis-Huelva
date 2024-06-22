pit <- res$cpo$pit
uniquant <- (1:138)/(138+1)
plot(uniquant, sort(pit), xlab="uniform quantiles", ylab="Sorted PIT values")
abline(0,1)

library(blockCV)

covariables_block_fold <- datos_huelva[,c("Distancia_Industria", "Distancia_Marsh", "Distancia_Urbano", 
                                          "Distancia_Agric", "Distancia_Phospho", "Distancia_Refineria", 
                                          "Distancia_Bare")]

datos_huelva_ras <- st_read("Datos/HuelvaResultados.shp")
zona <- raster::raster("Raster/tercero.tif")

# Leer el raster
zona <- raster("Raster/tercero.tif")

crs(zona) <- CRS("+init=epsg:25829")

sac2 <- blockCV::cv_spatial_autocor(x = datos_huelva_ras, column = "Pb__ppm_")

plot(sac2$variograms[[1]])

sb1 <- cv_spatial(x = datos_huelva_ras,
                  r = zona,
                  k = 5, # number of folds
                  size = 4725, # size of the blocks in metres
                  selection = "random", # random blocks-to-fold
                  iteration = 50)

cv_plot(cv = sb1,
        r = zona,
        raster_colors = terrain.colors(10, alpha = 0.5),
        label_size = 4) 

cv_plot(cv = sb1, 
        x = covariables_block_fold)


mae <- vector()
medae_k <- vector()
medae_k_normal <- vector()
elpd_total <- vector()
medae_k_ori <- vector()

for(k in 1:5){
  
  index_val=sb1$folds_list[[k]][[2]]
  index_train=sb1$folds_list[[k]][[1]]
  
  
}

for(k in 1:5){
  # Define the index
  
  index_val= sb1$folds_list[[k]][[2]]
  index_train= sb1$folds_list[[k]][[1]]
  
  max.edge = diff(range(st_coordinates(datos_huelva[index_train,])[,1]))/(3*3)
  
  bound.outer = diff(range(st_coordinates(datos_huelva[index_train,])[,1]))/2
  
  huelva.bdry <-  as(andalucia_huelva, "Spatial") %>% 
    inla.sp2segment()
  
  mesh_huelva <- inla.mesh.2d(boundary = huelva.bdry,
                              max.edge = c(0.7,0.9)*max.edge,
                              offset=c(max.edge, bound.outer),
                              cutoff = max.edge/3,
                              crs = 4326)
  # Define the weights
  A.train <- inla.spde.make.A(mesh=mesh_huelva, loc=coo[index_train,])
  A.val <- inla.spde.make.A(mesh=mesh_huelva, loc=coo[index_val,])
  
  spde <- inla.spde2.matern(mesh = mesh_huelva, alpha = 2, constr = TRUE)
  indexs <- inla.spde.make.index("s", spde$n.spde)
  
  # Build the stack
  covariates_list=c("ind","urban","mar","agric","refinery","phospho")
  
  stack.train <- inla.stack(tag = "train", data = list(y = datos_huelva$Hg__ppb_[index_train])
                            , A = list(1, A.train),
                            effects = list(data.frame(b0 = 1, 
                                                      ind = datos_huelva$Distancia_Industria[index_train],
                                                      mar = datos_huelva$Distancia_Marsh[index_train],
                                                      agric = datos_huelva$Distancia_Agric[index_train],
                                                      refi = datos_huelva$Distancia_Refineria[index_train],
                                                      urban = datos_huelva$Distancia_Urbano[index_train],
                                                      phospho = datos_huelva$Distancia_Phospho[index_train],
                                                      bare = datos_huelva$Distancia_Bare[index_train]),
                                           
                                           s = indexs))
  
  # The tag allows selected statistics for the desired sample to be extracted later
  stack.val <- inla.stack(data = list(pred=NA), # Set it to NA
                          A = list(1, A.val), 
                          effects = list(data.frame(b0 = 1, 
                                                    ind = grid_huelva$Industria[index_val],
                                                    mar = grid_huelva$Marsh[index_val],
                                                    agric = grid_huelva$Agric[index_val],
                                                    refi = grid_huelva$Refineria[index_val],
                                                    urban = grid_huelva$Urbano[index_val],
                                                    phospho = grid_huelva$Phospho[index_val],
                                                    bare = grid_huelva$Bare[index_val]
                          ), s = indexs),
                          tag = "val")
  
  join.stack <- inla.stack(stack.train, stack.val)
  
  # Write the test formula ####
  
  formula_test <- y ~ -1 + b0 + ind + bare + f(s, model = spde)
  
  # Fit the model
  spde_res <- inla(formula_test, family = "normal", 
                   data = inla.stack.data(join.stack),
                   control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
                   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                          return.marginals.predictor = TRUE,
                                          config = TRUE),
                   control.fixed=list(mean = 0,prec = 0.01,
                                      mean.intercept = 0,prec.intercept = 0.01))
  end_time=Sys.time()
  
  fun2 <- function()
    return(
      b0 +
        ind * grid_huelva$Industria[index_val] +
        bare * grid_huelva$Bare[index_val] +
        A.val %*% s
    )
  
  n_samples <- 1000
  
  samples <- inla.posterior.sample(n_samples, spde_res)
  
  prediccion_posterior_2 <- inla.posterior.sample.eval(fun2, samples)
  
  matrix_predictor <- matrix(NA, nrow = n_samples, ncol = length(index_val))
  vector_precision <- vector()
  
  for (i in 1:n_samples) {
    matrix_predictor[i, ] <- samples[[i]]$latent[index_val]
    vector_precision[i] <- samples[[i]]$hyperpar[[1]]
  }
  
  true_values = log(datos_huelva$Hg__ppb_[index_val])
  
  
  residu <-  sweep(matrix_predictor, 2, true_values) 
  
  log_dens <- dnorm(residu, sd=sqrt(1/vector_precision), log = TRUE)
  
  lpd <- apply(log_dens, 2, function (col) {logSumExp(col) - log(length(col))})
  
  sum(lpd)
  
  # Extract the fitted values
  index_inla_train = inla.stack.index(join.stack,"train")$data
  index_inla_val = inla.stack.index(join.stack,"val")$data
  
  results.train=spde_res$summary.fitted$mean[index_inla_train]
  results.val=spde_res$summary.fitted$mean[index_inla_val]
  
  
  
  print(median(abs(results.val - log(datos_huelva$Hg__ppb_[index_val]))))
  
  medae_k[k] <- median(abs(results.val - log(datos_huelva$Hg__ppb_[index_val])))
  medae_k_ori[k] <- median(abs(exp(results.val) - datos_huelva$Hg__ppb_[index_val]))
  medae_k_normal[k] <- median(abs(results.val-datos_huelva$Hg__ppb_[index_val]))
  
  print(1)
}

mean(medae_k)
mean(medae_k_ori)
mean(medae_k_normal)

originales <- c(47.64329, 40.81397, 108.6851, 14.77674, 75.36313)

sum(originales)
sqrt(138*var(originales))


fold_numbers <- 1:5 # Asumiendo que tienes 10 folds
data_for_plot <- data.frame(Fold = fold_numbers, MedAE = medae_k)


ggplot(data_for_plot, aes(x = Fold, y = MedAE)) +
  geom_line(group=1, colour="blue") +  # Línea que conecta los puntos de MedAE
  geom_point(aes(color = MedAE), size = 3) +  # Puntos de MedAE
  scale_color_gradient(low = "green", high = "red") +  # Colores para los puntos
  labs(title = "MedAE por Fold en la Validación Cruzada k-Fold",
       x = "Número de Fold",
       y = "Error Absoluto Mediano (MedAE)") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Oculta el título de la leyenda