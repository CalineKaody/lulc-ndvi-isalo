#############################################################################
# CHAINE LULC ISALO - VERSION 2 (contournement du bug sen2r / N9905)
#
# sen2r (abandonné depuis fin 2023) ne reconnaît pas le nouveau format de
# version de traitement ESA ("N9905") utilisé sur toutes les images
# actuelles. On contourne le problème en parlant DIRECTEMENT à Google Cloud
# via gsutil (déjà installé et authentifié chez toi via Google Cloud SDK),
# sans passer par sen2r pour la recherche/téléchargement.
#
# Les étapes 5, 6, 7 (classification, comparaison, stats) réutilisent les
# mêmes fonctions que le script original - rien à changer là-dessus.
#############################################################################

library(terra)
library(sf)
library(dplyr)

## ------------------------------------------------------------------------
## CONFIGURATION
## ------------------------------------------------------------------------
config <- list(
  aoi_path       = "data/ISALO.shp",  # chemin relatif : place le shapefile dans un dossier data/
  tuile          = "38KND",   # confirmé par l'erreur précédente
  periode_avant  = c("2019-01-01", "2019-06-30"),
  periode_apres  = c("2025-01-01", "2025-06-30"),
  n_classes      = 4,
  dossier_brut      = "sortie/01_brut",
  dossier_resultats = "sortie/05_resultats"
)
for (d in c(config$dossier_brut, config$dossier_resultats)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ------------------------------------------------------------------------
## 1. LISTER LES SCENES DISPONIBLES SUR LA TUILE (via gsutil)
## ------------------------------------------------------------------------
lister_scenes_tuile <- function(tuile) {
  utm   <- substr(tuile, 1, 2)
  bande <- substr(tuile, 3, 3)
  carre <- substr(tuile, 4, 5)
  chemin_bucket <- sprintf("gs://gcp-public-data-sentinel-2/L2/tiles/%s/%s/%s/",
                            utm, bande, carre)

  message("Interrogation de : ", chemin_bucket)
  resultat <- system2("gsutil", c("ls", chemin_bucket), stdout = TRUE, stderr = TRUE)

  # Ne garder que les lignes correspondant à des dossiers .SAFE/
  scenes <- grep("\\.SAFE/$", resultat, value = TRUE)
  if (length(scenes) == 0) {
    stop("Aucune scène trouvée. Vérifie que gsutil fonctionne (system2('gsutil', 'version'))",
         " et que la tuile '", tuile, "' est correcte.\nSortie brute :\n",
         paste(resultat, collapse = "\n"))
  }
  message(length(scenes), " scènes trouvées sur la tuile ", tuile)
  scenes
}

## ------------------------------------------------------------------------
## 2. FILTRER PAR PERIODE (date extraite du nom du fichier SAFE)
## ------------------------------------------------------------------------
filtrer_scenes_periode <- function(scenes, date_debut, date_fin) {
  noms <- basename(gsub("/$", "", scenes))
  # Format : S2X_MSIL2A_YYYYMMDDTHHMMSS_...
  dates_str <- sub("^S2[A-Z]_MSIL2A_([0-9]{8})T.*", "\\1", noms)
  dates <- as.Date(dates_str, format = "%Y%m%d")

  garder <- dates >= as.Date(date_debut) & dates <= as.Date(date_fin)
  scenes_filtrees <- scenes[garder]
  dates_filtrees <- dates[garder]

  if (length(scenes_filtrees) == 0) {
    stop("Aucune scène entre ", date_debut, " et ", date_fin,
         ". Essaie d'élargir la période dans config.")
  }

  # Trie par date, la plus récente en dernier
  ordre <- order(dates_filtrees)
  data.frame(
    scene = scenes_filtrees[ordre],
    date = dates_filtrees[ordre],
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------------------
## 3. TELECHARGER LES BANDES NECESSAIRES D'UNE SCENE (B04, B08, SCL)
## ------------------------------------------------------------------------
telecharger_bandes_scene <- function(chemin_scene, dossier_sortie) {
  dir.create(dossier_sortie, recursive = TRUE, showWarnings = FALSE)

  bandes_a_recuperer <- list(
    B04 = paste0(chemin_scene, "GRANULE/*/IMG_DATA/R10m/*_B04_10m.jp2"),
    B08 = paste0(chemin_scene, "GRANULE/*/IMG_DATA/R10m/*_B08_10m.jp2"),
    SCL = paste0(chemin_scene, "GRANULE/*/IMG_DATA/R20m/*_SCL_20m.jp2")
  )

  chemins_locaux <- list()
  for (b in names(bandes_a_recuperer)) {
    message("  Téléchargement bande ", b, "...")
    dest <- file.path(dossier_sortie, paste0(b, ".jp2"))
    res <- system2("gsutil", c("cp", shQuote(bandes_a_recuperer[[b]]), shQuote(dest)),
                    stdout = TRUE, stderr = TRUE)
    if (!file.exists(dest)) {
      warning("Echec du telechargement de la bande ", b, " :\n", paste(res, collapse = "\n"))
    } else {
      chemins_locaux[[b]] <- dest
    }
  }
  chemins_locaux
}

## ------------------------------------------------------------------------
## 4. PRETRAITER : NDVI + RECADRAGE SUR L'AOI + MASQUAGE NUAGES (SCL)
## ------------------------------------------------------------------------
pretraiter_scene <- function(chemins_bandes, aoi_sf) {
  b04 <- rast(chemins_bandes$B04)  # Rouge, 10m
  b08 <- rast(chemins_bandes$B08)  # PIR, 10m
  scl <- rast(chemins_bandes$SCL)  # Classification scène, 20m
  scl <- resample(scl, b04, method = "near")  # ramène le SCL à 10m

  ndvi <- (b08 - b04) / (b08 + b04)
  names(ndvi) <- "NDVI"

  # Masque nuages/ombres/cirrus/neige (classes SCL 3,8,9,10,11)
  classes_a_exclure <- c(3, 8, 9, 10, 11)
  masque_nuages <- !(scl %in% classes_a_exclure)
  ndvi_masque <- mask(ndvi, masque_nuages, maskvalue = 0)

  # Recadrage sur l'AOI Isalo
  aoi_proj <- st_transform(aoi_sf, crs(ndvi_masque))
  aoi_vect <- vect(aoi_proj)
  ndvi_final <- mask(crop(ndvi_masque, aoi_vect), aoi_vect)

  ndvi_final
}

## ------------------------------------------------------------------------
## 5-7. CLASSIFICATION, COMPARAISON, STATS (identiques au script original)
## ------------------------------------------------------------------------
classifier_occupation_sol <- function(ndvi_rast, n_classes = 4) {
  set.seed(184)
  valeurs <- values(ndvi_rast)
  complet <- complete.cases(valeurs)
  km <- kmeans(scale(valeurs[complet, , drop = FALSE]), centers = n_classes, nstart = 10)

  classif <- rast(ndvi_rast, nlyr = 1)
  vals_classif <- rep(NA, ncell(ndvi_rast))
  vals_classif[complet] <- km$cluster
  values(classif) <- vals_classif

  ndvi_moyen <- tapply(valeurs[complet, 1], km$cluster, mean, na.rm = TRUE)
  ordre <- order(ndvi_moyen)
  etiquettes <- c("Sol nu / Bati", "Vegetation clairsemee",
                   "Vegetation moderee", "Vegetation dense")[seq_len(n_classes)]
  correspondance <- setNames(etiquettes, names(ndvi_moyen)[ordre])

  list(raster = classif, correspondance = correspondance)
}

comparer_dans_le_temps <- function(classif_avant, classif_apres) {
  empile <- c(classif_avant, classif_apres)
  names(empile) <- c("avant", "apres")
  freq_tab <- crosstab(empile, useNA = FALSE)
  matrice <- as.data.frame(freq_tab)
  names(matrice) <- c("classe_avant", "classe_apres", "nb_pixels")
  carte_changement <- classif_avant != classif_apres
  list(matrice = matrice, carte = carte_changement)
}

produire_stats_surface <- function(classif_raster, correspondance, resolution_m = 10) {
  surface_par_pixel_ha <- (resolution_m^2) / 10000
  tab <- as.data.frame(freq(classif_raster))
  tab$classe <- correspondance[as.character(tab$value)]
  tab %>%
    group_by(classe) %>%
    summarise(surface_ha = sum(count) * surface_par_pixel_ha, .groups = "drop") %>%
    arrange(desc(surface_ha))
}

## ------------------------------------------------------------------------
## LANCEMENT ETAPE PAR ETAPE (exécute chaque bloc un par un dans la Console)
## ------------------------------------------------------------------------

# --- Etape A : charger l'AOI ---
aoi_sf <- st_read(config$aoi_path, quiet = TRUE)

# --- Etape B : lister toutes les scènes de la tuile (une seule fois) ---
toutes_scenes <- lister_scenes_tuile(config$tuile)

# --- Etape C : filtrer par période et choisir une scène par période ---
scenes_avant <- filtrer_scenes_periode(toutes_scenes, config$periode_avant[1], config$periode_avant[2])
scenes_apres <- filtrer_scenes_periode(toutes_scenes, config$periode_apres[1], config$periode_apres[2])
print(scenes_avant)   # regarde les dates disponibles et choisis-en une
print(scenes_apres)

# --- Etape D : télécharger la scène choisie pour chaque période ---
# (ici on prend la première par défaut - change l'indice [1] si tu préfères une autre date)
bandes_avant <- telecharger_bandes_scene(scenes_avant$scene[1], file.path(config$dossier_brut, "avant"))
bandes_apres <- telecharger_bandes_scene(scenes_apres$scene[1], file.path(config$dossier_brut, "apres"))

# --- Etape E : prétraiter (NDVI + masquage nuages + recadrage AOI) ---
ndvi_avant <- pretraiter_scene(bandes_avant, aoi_sf)
ndvi_apres <- pretraiter_scene(bandes_apres, aoi_sf)

plot(ndvi_avant, main = "NDVI - avant")
plot(ndvi_apres, main = "NDVI - apres")

# --- Etape F : classification ---
classif_avant <- classifier_occupation_sol(ndvi_avant, config$n_classes)
classif_apres <- classifier_occupation_sol(ndvi_apres, config$n_classes)

plot(classif_avant$raster, main = "Classification - avant")
plot(classif_apres$raster, main = "Classification - apres")

# --- Etape G : comparaison temporelle ---
changement <- comparer_dans_le_temps(classif_avant$raster, classif_apres$raster)
plot(changement$carte, main = "Zones de changement", col = c("grey90", "red"))

# --- Etape H : statistiques de surface ---
stats_avant <- produire_stats_surface(classif_avant$raster, classif_avant$correspondance)
stats_apres <- produire_stats_surface(classif_apres$raster, classif_apres$correspondance)
print(stats_avant)
print(stats_apres)

# --- Etape I : sauvegarder les résultats ---
writeRaster(classif_avant$raster, file.path(config$dossier_resultats, "classif_avant.tif"), overwrite = TRUE)
writeRaster(classif_apres$raster, file.path(config$dossier_resultats, "classif_apres.tif"), overwrite = TRUE)
write.csv(stats_avant, file.path(config$dossier_resultats, "stats_avant.csv"), row.names = FALSE)
write.csv(stats_apres, file.path(config$dossier_resultats, "stats_apres.csv"), row.names = FALSE)
message("Termine. Resultats dans : ", config$dossier_resultats)
