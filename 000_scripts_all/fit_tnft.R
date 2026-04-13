# Packages .............................................................................................................
.libPaths(c("/projappl/project_2004932/project_rpackages_4.4.0", .libPaths()))
if(!require("pacman")) {install.packages("pacman")}

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

pacman::p_load(
               tidyverse,
               datawizard,
               metaSEM)
 
i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))           

data <- read_csv("model_data_final.csv")

source("functions.R")

vars <- c(
  "TD_nestedness",
  "FD_turnover",
  "TP",
  "TDxTPmin",
  "Urban", "Agriculture", "Forest",
  "AT", "AP",
  "spatial_distance"
)

model <- '
  TP ~ Urban + Agriculture + Forest + AP

  TD_nestedness ~ TP + Urban + Agriculture + Forest +
                AT + AP + spatial_distance

  FD_turnover ~ TD_nestedness + TDxTPmin + TP +
                  Urban + Agriculture + Forest +
                  AT + AP + spatial_distance
                  
  TDxTPmin~~1*TDxTPmin
  Urban ~~ 1*Urban
  Agriculture ~~ 1*Agriculture
  Forest ~~ 1*Forest
  AT ~~ 1*AT
  AP ~~ 1*AP
  spatial_distance ~~ 1*spatial_distance
'

RAM <- lavaan2RAM(model, obs.variables = vars)

data_by_rep <- lapply(
  1:100,
  function(i) {
    data |>
      dplyr::filter(rep == i) |>
      mutate(TPmin = pmin(TP1,TP2)) |> 
      dplyr::select(
        dataset, Var1, Var2, rep,
        TD_nestedness, FD_turnover,
        spatial.distance,TPmin,
        TP, AT, AP,
        Urban_1500, Crops_1500, Forest_1500,
        Gamma_Diversity, FDGamma, Body_size, Group, Ecosystem, Dispersal_type
      ) |>
      as.data.frame()
  }
)

mods_by_ds <- data_by_rep[[1]] |>   # any rep is fine for dataset-level constants
   select(dataset, Gamma_Diversity, FDGamma, Body_size, Group, Ecosystem, Dispersal_type) |>
  
  group_by(dataset) |>
  
  summarise(
    Gamma_Diversity_log = first(log10(Gamma_Diversity)),
    FDGamma             = first(log10(FDGamma)),
    Body_size_log       = first(log10(Body_size)),
    Group               = as.character(first(na.omit(Group))),
    Ecosystem           = as.character(first(na.omit(Ecosystem))),
    Dispersal_type      = as.character(first(na.omit(Dispersal_type))),
    .groups = "drop"
  ) |>
  
  mutate(
    Dispersal_type = ifelse(Dispersal_type == "Seeds","Passive",Dispersal_type)
  ) |>
  mutate(
    Gamma_Diversity_log = as.numeric(standardise(Gamma_Diversity_log)),
    FDGamma             = as.numeric(standardise(FDGamma)),
    Body_size_log       = as.numeric(standardise(Body_size_log)),
    
    Group = factor(Group, levels = c("Consumer","Producer")),
    Ecosystem = factor(Ecosystem, levels = c("Lake","River")),
    Dispersal_type = factor(Dispersal_type, levels = c("Active","Passive")),
    
    Group = as.numeric(Group) - 1,
    Ecosystem = as.numeric(Ecosystem) - 1,
    Dispersal_type = as.numeric(Dispersal_type) - 1
  ) |>
  
  select(-dataset) |> 
  as.data.frame()

  data_i <- data_by_rep[[i]]   # small object now

  data_1 <- data_i |>
    mutate(
      Urban        = Urban_1500,
      Agriculture  = Crops_1500,
      Forest       = Forest_1500,
    ) |> 
    select(
      dataset, Var1, Var2, rep,
      TD_nestedness, FD_turnover,
      spatial.distance,TPmin,
      TP, AT, AP,
      Urban, Agriculture, Forest
    ) |>
    rename(spatial_distance = spatial.distance) |>
    mutate(
      rep     = factor(rep),
      dataset = factor(dataset)
    ) |>
    filter(rep == i) |>
    mutate(
      # GLOBAL mean-centering (no log, no scaling)
      TD_nestedness_c = TD_nestedness - mean(TD_nestedness, na.rm = TRUE),
      TPmin_c   = TPmin - mean(TPmin, na.rm = TRUE),
      
      # interaction term
      TDxTPmin = TD_nestedness_c * TPmin_c
    ) |>
    ungroup() |> 
    as.data.frame()
  
  by_ds <- data_1 |>
    group_by(dataset) |>
    summarise(
      R = list({
        X <- as.matrix(pick(all_of(vars)))
        R <- suppressWarnings(cor(X, use = "pairwise.complete.obs",method="spearman"))
        diag(R) <- 1
        R
      }),
      n = n_distinct(c(Var1, Var2)),
      .groups = "drop"
    ) 

  osma_dat <- Cor2DataFrame(by_ds$R, by_ds$n)
  
  osma_dat$data <- cbind(osma_dat$data, mods_by_ds)
  
  ## 4) Create Ax/Sx moderator matrices (multiple moderators -> list)
  mod_names <- colnames(mods_by_ds)
  
  idx_row <- function(v) which(rownames(RAM$A) == v)
  idx_col <- function(v) which(colnames(RAM$A) == v)
  
  tp_row <- idx_row("TP")
  td_row <- idx_row("TD_nestedness")
  fd_row <- idx_row("FD_turnover")
  
  tp_col   <- idx_col("TP")
  td_col   <- idx_col("TD_nestedness")
  dist_col <- idx_col("spatial_distance")
  
  env_cols <- which(colnames(RAM$A) %in% c("Urban","Agriculture","Forest"))
  
   Ax_list <- setNames(vector("list", length(mod_names)), mod_names)

for (m in mod_names) {

  base_mat <- create.modMatrix(RAM = RAM, output = "A", mod = m)
  temp_mat <- base_mat
  temp_mat[,] <- "0"

  ## ======================================================
  ## ECOSYSTEM TYPE
  ## modifies:
  ##  - landcover -> TP
  ##  - TP -> TD
  ##  - TP -> FD
  ## ======================================================
  if (grepl("Ecosystem", m)) {

    ## land cover -> TP
    temp_mat[tp_row, env_cols] <- base_mat[tp_row, env_cols]

    ## TP -> biodiversity
    temp_mat[td_row, tp_col] <- base_mat[td_row, tp_col]
    temp_mat[fd_row, tp_col] <- base_mat[fd_row, tp_col]
  }

  ## ======================================================
  ## DISPERSAL STRATEGY
  ## modifies:
  ##  - TP -> biodiversity
  ## ======================================================
  if (grepl("^Dispersal_type", m)) {

    temp_mat[td_row, tp_col] <- base_mat[td_row, tp_col]
    temp_mat[fd_row, tp_col] <- base_mat[fd_row, tp_col]
  }

  ## ======================================================
  ## TROPHIC GROUP
  ## modifies:
  ##  - TP -> biodiversity
  ## ======================================================
  if (grepl("Group", m)) {

    temp_mat[td_row, tp_col] <- base_mat[td_row, tp_col]
    temp_mat[fd_row, tp_col] <- base_mat[fd_row, tp_col]
  }

  ## ======================================================
  ## TAXONOMIC GAMMA DIVERSITY
  ## modifies:
  ##  - TP -> biodiversity
  ## ======================================================
  if (m == "Gamma_Diversity_log") {

    temp_mat[td_row, tp_col] <- base_mat[td_row, tp_col]
  }

  ## ======================================================
  ## FUNCTIONAL GAMMA
  ## modifies:
  ##  - TP -> biodiversity
  ## ======================================================
  if (m == "FDGamma") {
    temp_mat[fd_row, tp_col] <- base_mat[fd_row, tp_col]
  }

  ## ======================================================
  ## BODY SIZE
  ## modifies:
  ##  - TP -> biodiversity
  ## ======================================================
  if (m == "Body_size_log") {

    temp_mat[td_row, tp_col] <- base_mat[td_row, tp_col]
    temp_mat[fd_row, tp_col] <- base_mat[fd_row, tp_col]
  }

  Ax_list[[m]] <- temp_mat
}

  fit <- osmasem(
    model.name = "SN_to_FT_spearman",
    RAM  = RAM,
    Ax   = Ax_list,
    data = osma_dat,
    method = "REM",
  )
  

## Result: one row per estimated parameter per rep
## `coef_by_rep` contains ALL coefficients returned by coef(osmasem_fit),
## plus a `rep` column so you can summarise/plot later.

out_file <- glue::glue(
  "SNtoFT/metaSEM_rep{sprintf('%03d', i)}.rds"
)

saveRDS(fit, out_file)
