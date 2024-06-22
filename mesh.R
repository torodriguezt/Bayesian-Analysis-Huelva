library(dplyr)
library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(gridExtra)
library(INLA)
library(geoR)
library(fmesher)
library(lattice)
library(gridExtra)
library(readxl)
library(inlabru)
library(INLA)
library(splancs)

andalucia_huelva <- st_read("C:/Users/Tomas/Documents/INLA_Part2/Malla_no_water_2/malla.shp")
andalucia_huelva <- andalucia_huelva[,27]

andalucia_huelva %>% 
  st_set_crs(4326)


datos_huelva <- st_read("Datos/HuelvaResultados.shp")
datos_huelva <- st_zm(datos_huelva)
datos_huelva <- st_transform(datos_huelva, crs = 4326)

max.edge = diff(range(st_coordinates(datos_huelva)[,1]))/(3*3)

bound.outer = diff(range(st_coordinates(datos_huelva)[,1]))/3

huelva.bdry <-  as(andalucia_huelva, "Spatial") %>% 
  inla.sp2segment()

mesh_huelva <- inla.mesh.2d(boundary = huelva.bdry,
                            max.edge = c(0.7,0.9)*max.edge,
                            offset=c(max.edge, bound.outer),
                            cutoff = max.edge/3,
                            crs = 4326)

plot(mesh_huelva)


