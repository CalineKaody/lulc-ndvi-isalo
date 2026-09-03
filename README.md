# 🌍 Détection automatisée du changement d'occupation du sol — Isalo, Madagascar

Chaîne de traitement automatisée en **R** pour détecter les changements d'occupation du sol (LULC — Land Use/Land Cover) entre 2019 et 2025 sur la commune d'Isalo, à partir d'images satellite **Sentinel-2** (NDVI + classification non supervisée).

Projet réalisé dans le cadre du **M2 SIGD (Systèmes d'Information Géomatique et Décisionnel)** — EMIT.

## 🎯 Objectif

Automatiser un pipeline complet — de la recherche d'images satellite à la production de statistiques de surface — pour quantifier l'évolution de la végétation et de l'artificialisation des sols sur la commune d'Isalo entre deux périodes de référence.

## ⚙️ Méthodologie

1. **Recherche des scènes Sentinel-2** disponibles sur la tuile couvrant la zone, directement via Google Cloud Storage (bucket public `gcp-public-data-sentinel-2`)
2. **Filtrage temporel** : sélection d'une scène par période de référence (janvier–juin 2019 vs janvier–juin 2025)
3. **Téléchargement** des bandes nécessaires (B04 rouge, B08 proche infrarouge, SCL classification de scène) en 10m/20m de résolution
4. **Calcul du NDVI** `(B08 - B04) / (B08 + B04)`, avec masquage des nuages/ombres/cirrus via la bande SCL, puis recadrage sur la zone d'intérêt (AOI)
5. **Classification non supervisée** (k-means, 4 classes) du NDVI en catégories d'occupation du sol : sol nu/bâti, végétation clairsemée, modérée, dense
6. **Comparaison temporelle** entre les deux classifications (matrice de transition + carte binaire de changement)
7. **Statistiques de surface** par classe, en hectares, exportées en CSV

### 🔧 Point technique notable

`sen2r`, la librairie R habituellement utilisée pour automatiser la recherche/téléchargement Sentinel-2, est abandonnée depuis fin 2023 et ne reconnaît pas le format de version de traitement ESA actuel (`N9905`). La chaîne contourne ce blocage en interrogeant **directement le bucket Google Cloud public** via `gsutil`, sans dépendre de `sen2r` pour cette étape — le reste du pipeline (classification, comparaison, statistiques) reste inchangé.

## 📊 Résultats

| Classe | Surface 2019 (ha) | Surface 2025 (ha) |
|---|---:|---:|
| Sol nu / Bâti | 276 | 1 542 |
| Végétation clairsemée | 15 538 | 8 378 |
| Végétation modérée | 6 379 | 3 997 |
| Végétation dense | 5 923 | 3 302 |

La surface classée en **sol nu / bâti a été multipliée par plus de 5** entre les deux périodes, tandis que toutes les classes de végétation reculent. À noter : la couverture nuageuse diffère entre les deux scènes retenues (surface totale valide après masquage plus faible en 2025), ce qui doit être gardé à l'esprit dans l'interprétation des chiffres bruts — une comparaison en proportions relatives plutôt qu'en hectares absolus serait plus robuste pour une analyse plus poussée.

## 🖼️ Aperçu

- `captures/ndvi_avant.png` — carte NDVI, période "avant"
- `captures/classif_avant.png` / `classif_apres.png` — classifications par k-means
- `captures/changement.png` — carte binaire des zones de changement

## 🛠️ Technologies utilisées

| Composant | Outil |
|---|---|
| Langage | R |
| Traitement raster | `terra` |
| Données vectorielles | `sf` |
| Manipulation de données | `dplyr` |
| Accès aux images satellite | Google Cloud Storage (`gsutil`), bucket public Sentinel-2 |
| Classification | k-means (`stats`) |

## 🚀 Utilisation

### Prérequis
- R avec les packages `terra`, `sf`, `dplyr`
- [Google Cloud SDK](https://cloud.google.com/sdk) installé et `gsutil` accessible dans le PATH

### Étapes

```bash
git clone https://github.com/CalineKaody/lulc-ndvi-isalo.git
cd lulc-ndvi-isalo
```

Dans R :
```r
install.packages(c("terra", "sf", "dplyr"))
```

Le script `chaine_LULC_Isalo_v2.R` s'exécute étape par étape dans la console R (chaque bloc est commenté). Le fichier `data/ISALO.shp` (limite de la zone d'étude) est déjà inclus dans le dépôt.

## 📁 Structure du dépôt

```
lulc-ndvi-isalo/
├── chaine_LULC_Isalo_v2.R   # script principal
├── data/                     # shapefile de la zone d'étude (Isalo)
├── captures/                 # cartes NDVI, classifications, statistiques
└── rapport/                  # rapport complet du projet (.docx)
```

> Les rasters bruts (NDVI, différences) ne sont pas versionnés ici en raison de leur poids — disponibles sur demande.

## 👤 Auteure

FETISON Mioralalaina Câline — M2 SIGD, EMIT.
