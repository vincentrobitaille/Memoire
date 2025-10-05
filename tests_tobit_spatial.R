library(tidyr)
library(dplyr)
library(spatialreg)
library(spdep)

# jeu de données de test pour expérimenter
spData::afcon
data(lnd)
lnd

# Création de la matrice W (voisins) à partir des multipolygon
W <- lnd |> 
  poly2nb() |> 
  nb2mat()

