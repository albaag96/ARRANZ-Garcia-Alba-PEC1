if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")
if (!require("metabolomicsWorkbenchR", quietly = TRUE))
  install.packages("metabolomicsWorkbenchR")

library(SummarizedExperiment)
library(metabolomicsWorkbenchR)

# Descargar estudio de metabolomicsWorkbenchR
experiment = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST000002',
  output_item = 'SummarizedExperiment'
)

# Descripción de los datos
metadata(experiment)
colData(experiment)
colData(experiment)$Transplantation
rowData(experiment)
assay(experiment)




assay_data <- assay(experiment)
assay_data <- assay_data[, -1]


# Guardar el Summarized Experiment en formato binario
save(experiment, file = 'experiment.Rda')



if (!require("POMA", quietly = TRUE))
  BiocManager::install("POMA")

library(POMA)
library(ggtext)
library(magrittr)
library(readxl)
library(dplyr)
library(ggplot2)
library(readr)

# Modificar datos y cargar en POMA
metadata <- as.data.frame(colData(experiment))[, c(1, 6)]
metadata <- metadata %>%
  mutate(
    Transplantation = recode(Transplantation,
                         "Before transplantation" = "before",
                         "After transplantation" = "after")
  )
View(metadata)
View(assay(experiment))
data <- PomaCreateObject(metadata = metadata, features = t(assay(experiment))


# Detectar y eliminar valores perdidos
imputed <- data %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 20)

# Normalización
normalized <- imputed %>% 
  PomaNorm(method = "log_pareto")

# Visualización de datos no-normalizados y normalizados    

data$metada
PomaBoxplots(imputed, x = "samples")
PomaBoxplots(normalized, x = "samples")
PomaOutliers(normalized)$polygon_plot


# Análisis univariante
poma_uni <- PomaUnivariate(normalized)
poma_uni$result


# Análisis de componentes principales
poma_pca <- PomaPCA(pre_processed)
poma_pca$factors_plot
View(poma_pca$loadings)

