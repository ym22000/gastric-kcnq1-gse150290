Oui. Voici exactement ce qui a été fait.

## 1. Point de départ : le GSE

On est parti du jeu **`GSE150290`**.  
Dans le dossier gastric, j’ai téléchargé et utilisé les **matrices d’expression déjà traitées par échantillon** qui correspondent au compartiment **non immunitaire**.

Concrètement :
- les fichiers GSM traités sont dans `cache/`
- chaque fichier `.txt.gz` contient une **matrice gènes x cellules**
- les colonnes = cellules
- les lignes = gènes
- on a aussi gardé `raw_geo/` pour la traçabilité, mais la pipeline finale s’appuie surtout sur `cache/`

Donc ici, je n’ai pas refait l’alignement ou le comptage brut depuis fastq.  
Je suis parti des **matrices d’expression publiques déjà disponibles**.

## 2. Construction du t-SNE global

J’ai fusionné toutes les matrices non immunitaires en un seul objet Seurat.

On obtient :
- **13 022 cellules**
- **8 705 gènes**

Puis j’ai appliqué une pipeline Seurat proche de celle du papier :
- `NormalizeData`
- `FindVariableFeatures` avec `mean.var.plot`
- `mean.cutoff = c(0.0125, 6)`
- `dispersion.cutoff = c(0.5, Inf)`
- `ScaleData` avec régression de `nCount_RNA`
- `RunPCA(npcs = 30, seed.use = 12345)`
- `FindNeighbors(dims = 1:20)`
- `FindClusters(resolution = 0.8)`
- `RunTSNE(dims = 1:5, seed.use = 12345)`

## 3. Annotation des clusters globaux

Ensuite, les clusters globaux n’ont pas été annotés “à la main” juste en regardant le t-SNE.  
Je leur ai attribué un **programme de marqueurs moyen**.

J’ai défini plusieurs signatures simples :
- `tumor` : `EPCAM`, `CDH17`, `COL3A1`, `PDGFRB`
- `IM` : `TFF3`, `CDX1`, `CDX2`
- `PMC` : `GKN1`, `GKN2`, `MUC5AC`
- `GMC` : `MUC6`, `TFF2`
- `fibroblast` : `MMP2`, `PDGFRA`, `MYL9`, `FN1`, `CAV1`
- `EC` : `PLVAP`, `KDR`, `PTPRB`
- `enteroendocrine` : `CHGA`, `GAST`, `PROX1`

Pour chaque cluster global :
- j’ai calculé l’expression moyenne des gènes de chaque programme
- puis j’ai attribué au cluster le **lineage avec le score moyen le plus élevé**

C’est comme ça qu’on a obtenu les labels globaux du t-SNE :
- `GMC`
- `IM`
- `tumor`
- `fibroblast`
- `EC`
- `enteroendocrine`
- etc.

## 4. Comment les cellules tumorales ont été isolées

Après annotation du t-SNE global, j’ai pris **tous les clusters annotés `tumor`**.

Dans nos résultats, il y en avait **2** :
- cluster global `13`
- cluster global `3`

Ensemble, ils contenaient :
- **1408 cellules tumorales**

Donc ici, les cellules tumorales ont été isolées :
- **à partir du t-SNE global**
- **par annotation transcriptomique de cluster**
- et non par CNV ou inferCNV dans la version finale

C’est important :  
on n’a pas cherché à garder seulement les `400` cellules les plus conservatrices.  
On a préféré garder le **compartiment tumoral large**, pour mieux conserver l’hétérogénéité.

## 5. Construction du t-SNE tumoral

Une fois les **1408 cellules tumorales** extraites, j’ai recréé un nouvel objet Seurat avec **uniquement ces cellules**.

Puis j’ai refait exactement la même logique Seurat :
- `NormalizeData`
- `FindVariableFeatures(mean.var.plot, cutoff 0.0125–6, dispersion >= 0.5)`
- `ScaleData(vars.to.regress = nCount_RNA)`
- `RunPCA(npcs = 30, seed.use = 12345)`
- `FindNeighbors(dims = 1:20)`
- `FindClusters(resolution = 0.8)`
- `RunTSNE(dims = 1:5, seed.use = 12345)`

Ce t-SNE tumoral donne :
- **1408 cellules**
- **12 clusters tumoraux**

## 6. Comment on a défini intestinal vs diffuse dans le t-SNE tumoral

Ensuite, au lieu de forcer les clusters du papier, on a fait plus simple et plus robuste.

On a défini 2 pôles :
- **Intestinal** : `CDH17`, `REG4`, `MUC13`
- **Diffuse** : `IGFBP5`, `COL1A1`, `S100A4`, `TAGLN`, `EGR1`

Pour chaque cellule tumorale :
- j’ai calculé un **module score intestinal**
- et un **module score diffuse**

Puis :
- si `intestinal score >= diffuse score` → `Intestinal cells`
- sinon → `Diffuse cells`

C’est cette étape qui donne le t-SNE tumoral final annoté en :
- `Intestinal cells`
- `Diffuse cells`

## 7. Ce que ça donne au final

Dans le t-SNE tumoral :
- les clusters **1, 5, 6** sont surtout **intestinaux**
- les clusters **0, 2, 3, 4, 7, 8, 9, 10, 11** sont surtout **diffuse**

Et c’est sur cette base qu’on a ensuite projeté `KCNQ1` et testé sa corrélation.

## 8. Résumé ultra simple

En une phrase :

- on a pris le **GSE traité**
- construit un **t-SNE global non immunitaire**
- annoté les clusters par **signatures de marqueurs**
- extrait les **clusters globaux tumoraux**
- refait un **t-SNE tumoral**
- puis séparé les cellules tumorales en **pôle intestinal** versus **pôle diffuse**
