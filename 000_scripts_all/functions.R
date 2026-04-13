get_interaction_results <- function(models, facet = "tax") {
  
  ## ---- 1. TP effect (direção e magnitude média) ----
  tp_support <- get_sign_support(models) |>
    filter(term == "TP_log") |>
    transmute(
      tp_median   = median,
      tp_q2.5     = q2.5,
      tp_q97.5    = q97.5,
      tp_prop_pos = prop_pos
    )
  
    grids <- list(
      Urban_log  = quantile(model_data$Urban_log,  probs = seq(0.05, 0.95, length.out = 50), na.rm = TRUE),
      Agriculture_log  = quantile(model_data$Agriculture_log,  probs = seq(0.05, 0.95, length.out = 50), na.rm = TRUE),
      Forest_log = quantile(model_data$Forest_log, probs = seq(0.05, 0.95, length.out = 50), na.rm = TRUE)
    ) 
  
  
  
  ## ---- 2. Slopes ao longo dos gradientes ----
  slopes_all <- bind_rows(lapply(seq_along(models), function(i) {
    
    X <- models[[i]]
    
    rbind(
      # Urban
      as.data.frame(slopes(
        X,
        variables = "TP_log",
        by = "Urban_log",
        newdata = datagrid(Urban_log = grids[[1]]),
        vcov = FALSE,
        re.form = NA
      )) |>
        transmute(
          rep = i,
          moderator = "Urban_log",
          gradient = Urban_log,
          estimate = estimate
        ),
      
      # Crops
      as.data.frame(slopes(
        X,
        variables = "TP_log",
        by = "Agriculture_log",
        newdata = datagrid(Agriculture_log = grids[[2]]),
        vcov = FALSE,
        re.form = NA
      )) |>
        transmute(
          rep = i,
          moderator = "Agriculture_log",
          gradient = Agriculture_log,
          estimate = estimate
        ),
      
      # Forest
      as.data.frame(slopes(
        X,
        variables = "TP_log",
        by = "Forest_log",
        newdata = datagrid(Forest_log = grids[[3]]),
        vcov = FALSE,
        re.form = NA
      )) |>
        transmute(
          rep = i,
          moderator = "Forest_log",
          gradient = Forest_log,
          estimate = estimate
        )
    )
  }))
  
  ## ---- 3. Modulação do efeito de TP (min–max ao longo do gradiente) ----
  delta_slopes <- slopes_all |>
    group_by(rep, moderator) |>
    summarise(
      slope_min  = estimate[which.min(gradient)],
      slope_max  = estimate[which.max(gradient)],
      mean_slope = mean(estimate),
      delta_abs  = slope_max - slope_min,
      delta_pct  =  (slope_max - slope_min) / slope_min * 100,
      .groups = "drop"
    )
  
  ## ---- 4. Resumo final por moderador ----
  final_summary <- slopes_all |>
    
    # 1. Estatísticas por rep ao longo do gradiente
    group_by(rep, moderator) |>
    summarise(
      slope_min_rep = estimate[which.min(gradient)],
      slope_med_rep = estimate[which.min(abs(gradient - median(gradient)))],
      slope_max_rep = estimate[which.max(gradient)],
      
      # delta direcional ao longo do gradiente
      delta_pct_rep = (slope_max_rep - slope_min_rep) / slope_min_rep * 100,
      .groups = "drop"
    ) |>
    
    # 2. Agregar entre reps
    group_by(moderator) |>
    summarise(
      mean_slope_min    = mean(slope_min_rep, na.rm = TRUE),
      mean_slope_median = mean(slope_med_rep, na.rm = TRUE),
      mean_slope_max    = mean(slope_max_rep, na.rm = TRUE),
      
      median_delta_pct  = median(delta_pct_rep, na.rm = TRUE),
      
      prop_increase = mean(delta_pct_rep > 0, na.rm = TRUE),
      prop_decrease = mean(delta_pct_rep < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

get_sign_support <- function(m) {
  
  tidy_all <- lapply(m, tidy)
  
  tidy_fixef <- subset(
    bind_rows(tidy_all, .id = "rep"),
    effect == "fixed"
  )
  
  summary <- tidy_fixef %>%
    group_by(term) %>%
    dplyr::summarise(
      median   = median(estimate, na.rm = TRUE),
      q2.5     = quantile(estimate, 0.025, na.rm = TRUE),
      q97.5    = quantile(estimate, 0.975, na.rm = TRUE),
      prop_pos = mean(estimate > 0, na.rm = TRUE),
      prop_neg = mean(estimate < 0, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    mutate(direction =  ifelse(
      q2.5 > 0, "positive",
      ifelse(q97.5 < 0, "negative", "inconclusive")
    ))
  return(summary)
}

get_tp_slopes <- function(models, tp_values_raw, data) {

  tpmin_mean <- attr(data$TP_min_log, "center")
  tpmin_sd   <- attr(data$TP_min_log, "scale")

  tpmin_scaled <- (log10(tp_values_raw + 1) - tpmin_mean) / tpmin_sd
  
  slopes_all <-  purrr::map_dfr(seq_along(models), function(i) {
    
    marginaleffects::slopes(
      models[[i]],
      variables = "TP_log",
      newdata = datagrid(TP_min_log = tpmin_scaled),
      re.form = NA
    ) |>
      as.data.frame() |>
      dplyr::mutate(
        rep = i,
        TP_min_raw = tp_values_raw
      )
  })
}



# --- General wrapper to run TP_min interaction models by trophic_state ---

run_trophic_model <- function(data,
                              dimension = c("taxonomic", "functional"),
                              component = c("turnover", "nestedness")) {
  
  dimension <- match.arg(dimension)
  component <- match.arg(component)
#  data<-lake_within_band
  # Select response and extra predictor
  if (dimension == "taxonomic" && component == "turnover") {
    response <- "TD_turnover_log"
    extra    <- "Gamma_Diversity_log"
  }
  
  if (dimension == "functional" && component == "turnover") {
    response <- "FD_turnover_log"
    extra    <- "FDGamma"
  }
  
  if (dimension == "taxonomic" && component == "nestedness") {
    response <- "TD_nestedness_log"
    extra    <- "Gamma_Diversity_log"
  }
  
  if (dimension == "functional" && component == "nestedness") {
    response <- "FD_nestedness_log"
    extra    <- "FDGamma"
  }
  
  split_data <- dplyr::group_split(data, rep)
  
  lapply(split_data, function(x) {
    lmer(
      reformulate(
        termlabels = c(
          "spatial.distance",
          "TP_log * TP_min_log",
          "TP_log * (Urban_log + Agriculture_log + Forest_log)",
          "AT_log", "AP_log",
          "Group", 
          "Body_size_log", "Dispersal_type",
          extra, "(1|dataset)"
        ),
        response = response
      ) ,
      data = x,
      REML = FALSE
    )
  })
}


compute_sem_effects <- function(coef_vector) {
  # coef_vector: named numeric vector of SEM coefficients for ONE bootstrap
  # names like "AONB" meaning A -> B
  # ---- 1. Parse paths ----
  coef_names <- names(coef_vector)
  
  split_paths <- strsplit(coef_names, "ON")
  
  from <- sapply(split_paths, `[`, 2)
  to   <- sapply(split_paths, `[`, 1)
  
  paths <- data.frame(
    from = from,
    to   = to,
    beta = as.numeric(coef_vector),
    stringsAsFactors = FALSE
  )|>    dplyr::filter(!grepl("\\d", from))
  
  # ---- 2. Identify direct effects ----
  direct <- paths 
  
  # ---- 3. Identify indirect effects (2-step only) ----
  indirect_list <- list()
  k <- 1
  
  for (i in seq_len(nrow(paths))) {
    for (j in seq_len(nrow(paths))) {
      if (paths$to[i] == paths$from[j]) {
        indirect_list[[k]] <- data.frame(
          from = paths$from[i],
          to   = paths$to[j],
          beta = paths$beta[i] * paths$beta[j],
          via  = paths$to[i],
          stringsAsFactors = FALSE
        )
        k <- k + 1
      }
    }
  }
  
  indirect <- if (length(indirect_list) > 0) {
    do.call(rbind, indirect_list)
  } else {
    data.frame()
  }
  
  # ---- 4. Aggregate indirect effects ----
  if (nrow(indirect) > 0) {
    indirect_sum <- aggregate(
      beta ~ from + to,
      data = indirect,
      sum
    )
  } else {
    indirect_sum <- data.frame()
  }
  
  # ---- 5. Aggregate direct effects ----
  direct_sum <- aggregate(
    beta ~ from + to,
    data = direct,
    sum
  )
  
  # ---- 6. Compute total effects ----
  total <- merge(
    direct_sum,
    indirect_sum,
    by = c("from", "to"),
    all = TRUE,
    suffixes = c("_DE", "_IE")
  )
  
  total[is.na(total)] <- 0
  total$beta_TE <- total$beta_DE + total$beta_IE
  
  # ---- 7. Return tidy result ----
  list(
#    direct   = direct_sum,
#    indirect = indirect_sum,
    total    = total
  )
}

prepare_data_heatmap <- function(
    m,
    predictor_order = c(
      "Agriculture",
      "Urban",
      "Forest",
      "Precip",
      "Temp",
      "Spatial distance",
      "Total phosphorus",
      "Taxonomic turnover",
      "Taxonomic nestedness",
      "Taxonomic turnover × TPmin",
      "Taxonomic nestedness × TPmin"
    )
) {
  
  effects_decomposed <-
    dplyr::bind_rows(
      lapply(m, function(x)
        compute_sem_effects(coef(x))$total),
      .id = "boot"
    )
  
  effects_clean <- effects_decomposed %>%
    dplyr::mutate(
      from = dplyr::case_when(
        
        from == "Agriculture"      ~ "Agriculture",
        from == "Urban"            ~ "Urban",
        from == "Forest"           ~ "Forest",
        from == "AP"               ~ "Precip",
        from == "AT"               ~ "Temp",
        from == "TP"               ~ "Total phosphorus",
        from == "TD_turnover"      ~ "Taxonomic turnover",
        from == "TD_nestedness"    ~ "Taxonomic nestedness",
        from == "FD_turnover"      ~ "Functional turnover",
        from == "FD_nestedness"    ~ "Functional nestedness",
        from == "spatial_distance" ~ "Spatial distance",
        
        ## ambiguous interaction resolved from the model content
        from == "TDxTPmin" & any(.data$from == "TD_turnover")   ~
          "Taxonomic turnover × TPmin",
        
        from == "TDxTPmin" & any(.data$from == "TD_nestedness") ~
          "Taxonomic nestedness × TPmin",
        
        TRUE ~ from
      ),
      
      to = dplyr::recode(
        to,
        "TD_turnover"   = "Taxonomic turnover",
        "TD_nestedness" = "Taxonomic nestedness",
        "FD_turnover"   = "Functional turnover",
        "FD_nestedness" = "Functional nestedness",
        "TP"            = "Total phosphorus"
      )
      
    ) |> 
    dplyr::mutate(
      from = factor(from, levels = predictor_order)
    )
  
  osmasem_standardised_effects <- effects_clean |>
    tidyr::pivot_longer(
      !c(from, to, boot),
      names_to  = "effect_type",
      values_to = "estimate"
    ) |>
    dplyr::mutate(
      effect_type = gsub("beta_", "", effect_type),
      effect_type = dplyr::recode(
        effect_type,
        "DE" = "Direct effect",
        "IE" = "Indirect effect",
        "TE" = "Total effect"
      )
    ) |>
    dplyr::group_by(from, to, effect_type) |>
    dplyr::summarise(
      median  = mean(estimate, na.rm = TRUE),
      support = 100 * mean(
        sign(estimate) == sign(mean(estimate, na.rm = TRUE)),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(effect_type != "Total effect") |>
    dplyr::mutate(
      effect_type = forcats::fct_relevel(
        effect_type,
        "Direct effect",
        "Indirect effect"
      ),
      to = forcats::fct_relevel(
        to,
        "Total phosphorus",
        "Taxonomic turnover",
        "Taxonomic nestedness",
        "Functional turnover",
        "Functional nestedness"
      )
      )
  
  osmasem_standardised_effects
}


plot_osmasem_path <- function(osmasem_df,
                              support_min = 90,
                              small_cut = 0.02) {
  
  stopifnot(
    requireNamespace("dplyr"),
    requireNamespace("ggplot2"),
    requireNamespace("ggforce"),
    requireNamespace("funkyheatmap"),
    requireNamespace("colorspace"),
    requireNamespace("tibble"),
    requireNamespace("grid")
  )
  
  ## ------------------------------------------------------------
  ## Fixed nodes
  ## ------------------------------------------------------------

  nodes <- tribble(
    ~id,        ~x, ~y, ~shape, ~w,  ~h,  ~fill,     ~label,
    "Land",     1.8, 10, "rect",  1.6, 0.9, "#D9D9E8", "Land Cover",
    "Climate",  4.5, 10, "rect",  1.6, 0.9, "#D9D9E8", "Climate",
    "TP",       3, 6, "round", 2, 0.9, "#5F7282", "Total phosphorus",
    "Spatial",  3.125, 1, "rect", 1, 1.0, "#BFE3A8", "Spatial\ndistance",
    "Tax",      1.25, 3.0, "hex", 1.4, NA, "#92C5DEFF", "Taxonomic\nturnover",
    "Func",     5, 3.0, "hex", 1.6, NA, "#F9AE78", "Functional\nnestedness"
  )
  

  # helpers
  rect_nodes  <- filter(nodes, shape == "rect")
  round_nodes <- filter(nodes, shape == "round")
  hex_nodes   <- filter(nodes, shape == "hex")
  
  
  
  rect_nodes  <- dplyr::filter(nodes, shape == "rect")
  round_nodes <- dplyr::filter(nodes, shape == "round")
  hex_nodes   <- dplyr::filter(nodes, shape == "hex")
  
  ## ------------------------------------------------------------
  ## Fixed bent edges (routing)
  ## ------------------------------------------------------------
  edges_bent <- tibble::tribble(
    ~edge, ~x,  ~y,
    
    "Land_TP", 2.2, 9.55,
    "Land_TP", 2.2, 6.66,
    
    "Land_Tax", 1.1, 9.55,
    "Land_Tax", 1.1, 3.75,
    
    "Land_Func", 2, 10,
    "Land_Func", 2, 10.6,
    "Land_Func", 6, 10.6,
    "Land_Func", 6, 3,
    "Land_Func", 5.8, 3,
    
    "Climate_TP", 3.8, 9.55,
    "Climate_TP", 3.8, 6.6,
    
    "Climate_Tax", 4, 10,
    "Climate_Tax", 4, 10.7,
    "Climate_Tax", 0, 10.7,
    "Climate_Tax", 0, 3,
    "Climate_Tax", .45, 3,
    
    "Climate_Func", 5, 10,
    "Climate_Func", 5, 3.7,
    
    "TP_Tax", 2.2, 6,
    "TP_Tax", 2.2, 4.5,
    "TP_Tax", 1.65, 3.7,
    
    "TP_Func", 3.8, 6,
    "TP_Func", 3.8, 4.5,
    "TP_Func", 4.59, 3.65,
    
    "Spatial_Tax", 2.7, 1,
    "Spatial_Tax", 1.2, 1,
    "Spatial_Tax", 1.2, 2.3,
    
    "Spatial_Func", 3.5, 1,
    "Spatial_Func", 4.8, 1,
    "Spatial_Func", 4.8, 2.3,
    
    "Tax_Func", 1.7, 3,
    "Tax_Func", 4.25, 3,
    
    "Tax_Func_int", 1.6, 3.3,
    "Tax_Func_int", 2.2, 4.2,
    "Tax_Func_int", 3.77, 4.2,
    "Tax_Func_int", 4.45, 3.44
  )
  
  ## ------------------------------------------------------------
  ## Semantic routing map
  ## ------------------------------------------------------------
  edge_key <- tibble::tribble(
    ~from, ~to_role, ~edge_id, ~label_x, ~label_y, ~group_id,
    
    "Agriculture","Total phosphorus","Land_TP",2.2,8.4,"Land_TP",
    "Urban","Total phosphorus","Land_TP",2.2,8.4,"Land_TP",
    "Forest","Total phosphorus","Land_TP",2.2,8.4,"Land_TP",
    
    "Agriculture","Functional","Land_Func",6.0,8.0,"Land_Func",
    "Urban","Functional","Land_Func",6.0,8.0,"Land_Func",
    "Forest","Functional","Land_Func",6.0,8.0,"Land_Func",
    
    "Precip","Total phosphorus","Climate_TP",4.0,8.1,"Climate_TP",
    "Temp","Total phosphorus","Climate_TP",4.0,8.1,"Climate_TP",
    
    "Agriculture","Taxonomic","Land_Tax",1.1,6.8,"Land_Tax",
    "Forest","Taxonomic","Land_Tax",1.1,6.8,"Land_Tax",
    
    "Precip","Taxonomic","Climate_Tax",0.0,5.3,"Climate_Tax",
    "Temp","Taxonomic","Climate_Tax",0.0,5.3,"Climate_Tax",
    
    "Precip","Functional","Climate_Func",5.0,5.7,"Climate_Func",
    "Temp","Functional","Climate_Func",5.0,5.7,"Climate_Func",
    
    "Total phosphorus","Taxonomic","TP_Tax",2.2,5.2,"TP_Tax",
    "Total phosphorus","Functional","TP_Func",3.8,5.2,"TP_Func",
    
    "Spatial distance","Taxonomic","Spatial_Tax",2.0,1.0,"Spatial",
    "Spatial distance","Functional","Spatial_Func",4.0,1.0,"Spatial",
    
    "Taxonomic turnover","Functional","Tax_Func",3.0,3.0,"Tax_Func",
    "Taxonomic nestedness","Functional","Tax_Func",3.0,3.0,"Tax_Func",
    
    "Taxonomic turnover × TPmin","Functional",
    "Tax_Func_int",3.0,4.2,"Tax_Func_int",
    
    "Taxonomic nestedness × TPmin","Functional",
    "Tax_Func_int",3.0,4.2,"Tax_Func_int"
  )
  
  ## ------------------------------------------------------------
  ## Only supported direct effects drive arrows
  ## ------------------------------------------------------------
  df <- osmasem_df |>
    dplyr::filter(effect_type == "Direct effect",
                  support >= support_min) |>
    dplyr::mutate(
      to_role = dplyr::case_when(
        grepl("^Taxonomic", to)  ~ "Taxonomic",
        grepl("^Functional", to) ~ "Functional",
        TRUE                     ~ to
      )
    ) |>
    dplyr::left_join(edge_key, by = c("from", "to_role")) |>
    dplyr::filter(!is.na(edge_id))
  
  tax_label  <- unique(df$to[grepl("^Taxonomic", df$to)])
  func_label <- unique(df$to[grepl("^Functional", df$to)])
  
  tax_label  <- tax_label[1]
  func_label <- func_label[1]
  
  nodes <- nodes |>
    dplyr::mutate(
      label = dplyr::case_when(
        id == "Tax"  ~ gsub(" ", "\n", tax_label),
        id == "Func" ~ gsub(" ", "\n", func_label),
        TRUE ~ label
      )
    )
  
  ## categorical blocks rule
  group_strength <- df |>
    dplyr::group_by(group_id) |>
    dplyr::summarise(max_abs = max(abs(median)), .groups="drop")
  
  df <- df |>
    dplyr::left_join(group_strength, by="group_id") |>
    dplyr::mutate(strong_group = max_abs >= small_cut)
  
  ## ------------------------------------------------------------
  ## colour logic
  ## ------------------------------------------------------------

  
  df <- df |>
    dplyr::mutate(
      
      arrow_col = dplyr::case_when(
        median >= 0 & abs(median) < small_cut  ~ colorspace::lighten("#0072B2", 0.6),
        median >= 0 & abs(median) >= small_cut ~ "#0072B2",
        
        median <  0 & abs(median) < small_cut  ~ colorspace::lighten("#C9483B", 0.6),
        median <  0 & abs(median) >= small_cut ~ "#C9483B"
      ),
      
      label_col = ifelse(abs(median) < small_cut, "grey50", "black")
    )
  
  ## one colour per routed edge
  edge_cols <- df |>
    dplyr::group_by(edge_id) |>
    dplyr::slice_max(abs(median), n=1, with_ties=FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(edge_id, arrow_col)
  
  edges_plot <- edges_bent |>
    dplyr::left_join(edge_cols,
                     by = c("edge"="edge_id")) |>
    dplyr::filter(!is.na(arrow_col))
  
  ## labels
  edge_labels <- df |>
    dplyr::group_by(edge_id, label_x, label_y, group_id) |>
    dplyr::reframe(
      
      label = {
        
        if (group_id[1] %in% c("Land_TP","Land_Tax","Land_Func",
                               "Climate_TP","Climate_Tax","Climate_Func")) {
          
          paste0(from, ": ", sprintf("%.2f", median), collapse = "\n")
          
        } else if (group_id[1] == "Tax_Func_int") {
          
          paste0("× TPmin\n", sprintf("%.2f", median[1]))
          
        } else {
          
          sprintf("%.2f", median[which.max(abs(median))])
          
        }
      },
      
      col = if (group_id[1] %in% c("Land_TP","Land_Tax","Land_Func",
                                   "Climate_TP","Climate_Tax","Climate_Func")) {
        ifelse(any(strong_group), "black", "grey50")
      } else {
        ifelse(strong_group[which.max(abs(median))], "black", "grey50")
      }
      
    ) |>
    dplyr::rename(x = label_x, y = label_y)
  
  edge_labels <- edge_labels |>
    dplyr::left_join(
      edge_cols,
      by = c("edge_id" = "edge_id")
    )
  
  ## ------------------------------------------------------------
  ## plot
  ## ------------------------------------------------------------
  ggplot2::ggplot() +
    
    ggplot2::geom_path(
      data = edges_plot,
      ggplot2::aes(x, y, group=edge, colour=arrow_col),
      linewidth=.6,
      lineend="round",
      arrow=ggplot2::arrow(type="closed",
                           length=grid::unit(3,"mm"))
    ) +
    
    ggplot2::geom_rect(
      data=rect_nodes,
      ggplot2::aes(xmin=x-w/2, xmax=x+w/2,
                   ymin=y-h/2, ymax=y+h/2),
      fill=rect_nodes$fill, 
      colour = colorspace::darken(rect_nodes$fill, 0.3)
    ) +
    
    funkyheatmap::geom_rounded_rect(
      data=round_nodes,
      ggplot2::aes(xmin=x-w/2, xmax=x+w/2,
                   ymin=y-h/2, ymax=y+h/2),
      radius=grid::unit(4,"mm"),
      fill=round_nodes$fill, 
      colour = colorspace::darken(round_nodes$fill, 0.3)
    ) +
    
    ggforce::geom_regon(
      data=hex_nodes,
      ggplot2::aes(x0=x,
                   y0=y,
                   r=.7,
                   sides=6,
                   angle=0,
                   fill=fill,
                   colour = colorspace::darken(fill, 0.3)),
      
    ) +
    
    ggplot2::scale_fill_identity() +
    
    ggplot2::geom_text(
      data=nodes,
      ggplot2::aes(x=x,y=y,label=label),
      colour=c(rep("black",2),"white","black","black","black"),
      size=4
    ) +
    
    ggplot2::geom_label(
      data = edge_labels,
      ggplot2::aes(
        x = x, y = y,
        label = label,
        colour = col,
      ),
      border.colour = edge_labels$arrow_col,
      fill = "white",
      label.size = .3,
      size = 3.6,
      label.r = grid::unit(0,"pt"),
      show.legend = FALSE
    )+
    
    ggplot2::scale_colour_identity() +
    ggplot2::coord_cartesian(xlim=c(-1,7), ylim=c(0,11), expand=FALSE) +
    ggplot2::theme_void()
}

extract_osmasem_moderators <- function(osma_obj){
  
  M <- osma_obj$Mmatrix
  ax_names <- grep("^Ax[0-9]+$", names(M), value = TRUE)
  
  mods <- lapply(ax_names, function(ax){
    labs <- M[[ax]]@labels
    labs <- unique(as.vector(labs))
    labs <- labs[!is.na(labs)]
    labs <- sub("^data\\.","",labs)
    if(length(labs)==0) return(NA_character_)
    labs[1]
  })
  
  idx <- as.numeric(sub("^Ax","",ax_names))
  ord <- order(idx)
  
  mods <- unlist(mods[ord],use.names=FALSE)
  names(mods) <- paste0("_",idx[ord])
  mods
}

extract_osmasem_paths <- function(model){
  
  coefs <- stats::coef(model)
  mods  <- extract_osmasem_moderators(model)
  
  tibble(
    term = names(coefs),
    estimate = as.numeric(coefs)
  ) |>
    
    filter(str_detect(term,"_[0-9]+$")) |>
    
    mutate(
      moderator_id = str_extract(term,"_[0-9]+$"),
      moderator    = unname(mods[moderator_id]),
      path         = str_remove(term,"_[0-9]+$"),
      response     = str_match(path,"^(.*?)ON(.*)$")[,2],
      predictor    = str_match(path,"^(.*?)ON(.*)$")[,3]
    ) |> 
    
    filter(!is.na(moderator)) |>
    
    mutate(
      
      resp_code = short_code(response),
      pred_code = short_code(predictor),
      
      moderator_label = dplyr::case_when(
        moderator == "Group" ~ "Consumer vs Producer",
        moderator == "Dispersal_type" ~ "Active vs Passive",
        moderator == "Ecosystem" ~ "Lake vs River",
        moderator == "Gamma_Diversity_log" ~ "Taxonomic γ-diversity",
        moderator == "FDGamma" ~ "Functional γ-diversity",
        moderator == "Body_size_log" ~ "Body size",
        TRUE ~ moderator
      ),
      
      path_label = case_when(
        
        predictor=="Urban" ~ "Urban → TP",
        predictor=="Agriculture" ~ "Agriculture → TP",
        predictor=="Forest" ~ "Forest → TP",
        
        predictor=="TP" ~ paste0("TP → ",nice_label(resp_code)),
        
        TRUE ~ paste0(resp_code,":",pred_code)
      )
    )
}

short_code <- function(x){
  dplyr::case_when(
    x == "TD_turnover"   ~ "TT",
    x == "TD_nestedness" ~ "TN",
    x == "FD_turnover"   ~ "FT",
    x == "FD_nestedness" ~ "FN",
    x == "TP" ~ "TP",
    x == "spatial_distance" ~ "SD",
    TRUE ~ x
  )
}

nice_label <- function(code){
  dplyr::case_when(
    code == "TT" ~ "Tax. turn.",
    code == "TN" ~ "Tax. nest.",
    code == "FT" ~ "Func. turn.",
    code == "FN" ~ "Func. nest.",
    TRUE ~ code
  )
}

plot_moderators <- function(m, conf = 0.95){

  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tibble)
    library(purrr)
    library(ggplot2)
    library(ggforce)
  })
  


  
  boot_df <- purrr::imap_dfr(
    m,
    ~extract_osmasem_paths(.x) |>
      mutate(draw=.y,.before=1)
  )
  
  alpha <- (1-conf)/2
  
  forest_df <- boot_df |>
    
    group_by(
      path_label,
      moderator,
      moderator_label
    ) |>
    
    summarise(
      median = mean(estimate),
      lwr    = quantile(estimate,alpha),
      upr    = quantile(estimate,1-alpha),
      .groups="drop"
    ) |>
    
    filter(str_detect(path_label,"TP"))
  
  ## facet order requested
  forest_df$moderator_label <- factor(
    forest_df$moderator_label,
    levels=c(
      "Lake vs River",
      "Active vs Passive",
      "Consumer vs Producer",
      "Body size",
      "Taxonomic γ-diversity",
      "Functional γ-diversity"
    )
  )
  
  ## ecosystem facet internal ordering
  eco_order <- c(
    "Forest → TP",
    "Agriculture → TP",
    "Urban → TP",
    "TP → Tax. turn.",
    "TP → Tax. nest.",
    "TP → Func. turn.",
    "TP → Func. nest."
  )
  
  forest_df$path_label <- factor(
    forest_df$path_label,
    levels=eco_order
  )
  
  ggplot(
    forest_df,
    aes(x=median,y=path_label)
  )+
    
    geom_vline(
      xintercept=0,
      linetype=2,
      colour="grey40"
    )+
    
    geom_errorbarh(
      aes(xmin=lwr,xmax=upr),
      height=0,
      linewidth=.4
    )+
    
    geom_point(size=2.6)+
    
    ggforce::facet_col(
      ~moderator_label,
      scales="free_y",
      space="free"
    )+
    
    labs(
      x="Moderation effect",
      y=NULL
    )+
    
    theme_bw(base_size=11)+
    
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill="grey90"),
      strip.text=element_text(face="bold"),
      panel.spacing.y=grid::unit(2,"mm"),
      plot.margin=margin(3,3,3,3),
      legend.position="none"
    )
}

predict_surface <- function(model_list, land_var, land_grid){

  pred_df <- purrr::imap_dfr(model_list, function(mod, rep_id){
    
    base <- marginaleffects::datagrid(model = mod)
    
    # create grid of TP × land use
    grid <- expand.grid(
      TP_log = tp_grid,
      gradient = land_grid
    )
    
    # replicate base row for each grid point
    nd <- base[rep(1, nrow(grid)), ]
    
    nd$TP_log <- grid$TP_log
    nd[[land_var]] <- grid$gradient
    
    marginaleffects::predictions(
      mod,
      newdata = nd,
      re.form = NA,
      transform = exp
    ) |>
      dplyr::select(estimate, TP_log, dplyr::all_of(land_var)) |>
      dplyr::mutate(rep = rep_id)
    
  })
  
  pred_df |>
    dplyr::group_by(
      TP_log,
      gradient = .data[[land_var]]
    ) |>
    dplyr::summarise(
      estimate = mean(estimate),
      .groups = "drop"
    )
}









extract_fd_effects_one_model <- function(eff_df, model_name = NULL) {
  
  req_cols <- c("from", "to", "effect_type", "median", "support")
  miss_cols <- setdiff(req_cols, names(eff_df))
  if (length(miss_cols) > 0) {
    stop("eff_df is missing required columns: ",
         paste(miss_cols, collapse = ", "))
  }
  dat <- eff_df %>%
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      effect_type = as.character(effect_type)
    )
  
  # 找功能响应变量
  fd_targets <- unique(dat$to[grepl("^Functional\\s", dat$to)])
  
  if (length(fd_targets) == 0) {
    stop("No functional target found in column `to`.")
  }
  
  if (length(fd_targets) > 1) {
    warning("Multiple functional targets found: ",
            paste(fd_targets, collapse = ", "),
            ". Using the first one: ", fd_targets[1])
  }
  
  fd_target <- fd_targets[1]
  
  out <- dat %>%
    dplyr::filter(to == fd_target) %>%
    dplyr::mutate(
      effect_type = dplyr::recode(
        effect_type,
        "Direct effect"   = "DE",
        "Indirect effect" = "IE"
      ),
      model = if (is.null(model_name)) fd_target else model_name,
      fd_target = fd_target
    ) %>%
    dplyr::select(from, effect_type, median, support, model, fd_target)
  
  return(out)
}






