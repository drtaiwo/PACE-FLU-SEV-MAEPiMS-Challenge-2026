# =============================================================================
# PACE-FLU FINAL M4
# DEFINITIVE LEAKAGE-FREE 24-ORIGIN BAYESIAN VALIDATION
#
# FINAL FROZEN MODEL:
#
#   Y_{i,t+1} ~ Negative Binomial(mu_{i,t+1}, phi)
#
#   log(mu_{i,t+1}) =
#
#       log(P_i / 100000)
#       + alpha
#       + alpha_i
#       + beta_M  * Memory_it
#       + beta_N  * National_t
#       + beta_B  * Borrowing_it
#       + beta_P  * Phase_it
#       + beta_BP * (Borrowing_it * Phase_it)
#
#
# IMPORTANT:
#
# 1. M4 architecture is frozen.
#
# 2. No intensity extension is added.
#
# 3. borrow_phase = borrowing * p_growth
#    and is NOT separately standardized.
#
# 4. Current-season training predictor rows stop at origin - 1.
#
# 5. The forecast row at week t predicts week t+1.
#
# 6. Previous seasons are fully available.
#
# 7. Connectivity and growth scale are learned using training data only.
#
# 8. Connectivity predictor rows from the current season therefore never
#    use the held-out t+1 outcome.
#
# 9. Numerical settings:
#       init = 0
#       step_size = 0.01
#       adapt_delta = 0.99
#       max_treedepth = 15
#
# 10. This is the FINAL M4 validation.
#
# =============================================================================



# =============================================================================
# 0. PACKAGES
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(cmdstanr)
library(posterior)
library(ggplot2)



# =============================================================================
# 1. GLOBAL SETTINGS
# =============================================================================

SEED <- 20260925

set.seed(SEED)


STATE_DATA_FILE <-
  "data/nigeria_flu_weekly_by_state.csv"


OUTPUT_DIR <-
  "PACE_FLU_FINAL_M4_24_ORIGIN_VALIDATION"


FIT_DIR <-
  file.path(
    OUTPUT_DIR,
    "fits"
  )


CHECKPOINT_DIR <-
  file.path(
    OUTPUT_DIR,
    "checkpoints"
  )


PLOT_DIR <-
  file.path(
    OUTPUT_DIR,
    "plots"
  )


STAN_DIR <-
  file.path(
    OUTPUT_DIR,
    "stan"
  )


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  FIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  CHECKPOINT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  PLOT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  STAN_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)



# =============================================================================
# 2. SEASONS AND ORIGINS
# =============================================================================

SEASON_LEVELS <- c(
  "2023/2024",
  "2024/2025",
  "2025/2026"
)


ORIGINS <- c(
  20L,
  24L,
  28L,
  32L,
  36L,
  40L,
  44L,
  48L
)


VALIDATION_GRID <-
  tidyr::expand_grid(
    
    season =
      SEASON_LEVELS,
    
    origin =
      ORIGINS
  ) %>%
  
  dplyr::mutate(
    
    season_index =
      match(
        .data$season,
        SEASON_LEVELS
      ),
    
    target_week =
      .data$origin + 1L
  ) %>%
  
  dplyr::arrange(
    .data$season_index,
    .data$origin
  )



# =============================================================================
# 3. MEMORY SETTINGS
# =============================================================================

L_MEMORY <- 8L


MEMORY_DECAY <- 0.35


MEMORY_WEIGHTS <-
  exp(
    -MEMORY_DECAY *
      (0:(L_MEMORY - 1L))
  )


MEMORY_WEIGHTS <-
  MEMORY_WEIGHTS /
  sum(
    MEMORY_WEIGHTS
  )



# =============================================================================
# 4. CONNECTIVITY SETTINGS
# =============================================================================

MIN_EDGE_N <- 20L



# =============================================================================
# 5. MCMC SETTINGS
# =============================================================================

N_CHAINS <- 4L


N_PARALLEL_CHAINS <- 4L


N_WARMUP <- 1000L


N_SAMPLING <- 1000L


ADAPT_DELTA <- 0.99


MAX_TREEDEPTH <- 15L


STAN_INIT <- 0


INITIAL_STEP_SIZE <- 0.01


MAX_PRED_DRAWS <- 4000L



# =============================================================================
# 6. GENERAL FUNCTIONS
# =============================================================================

safe_print <- function(x) {
  
  print(
    as.data.frame(x),
    row.names = FALSE
  )
}



normalize_season <- function(x) {
  
  x <- as.character(x)
  
  dplyr::case_when(
    
    x %in%
      c(
        "2023/24",
        "2023/2024"
      ) ~
      "2023/2024",
    
    x %in%
      c(
        "2024/25",
        "2024/2025"
      ) ~
      "2024/2025",
    
    x %in%
      c(
        "2025/26",
        "2025/2026"
      ) ~
      "2025/2026",
    
    TRUE ~
      NA_character_
  )
}



safe_sd <- function(x) {
  
  ans <-
    stats::sd(
      x,
      na.rm = TRUE
    )
  
  if (
    !is.finite(ans) ||
    ans <= 1e-8
  ) {
    
    ans <- 1
  }
  
  ans
}



check_finite <- function(
    x,
    name
) {
  
  if (
    any(
      !is.finite(x)
    )
  ) {
    
    stop(
      paste0(
        "Non-finite values detected in ",
        name
      )
    )
  }
}



safe_percent_improvement <- function(
    baseline,
    alternative
) {
  
  ifelse(
    
    is.finite(baseline) &
      abs(baseline) >
      1e-12,
    
    100 *
      (
        baseline -
          alternative
      ) /
      baseline,
    
    NA_real_
  )
}



# =============================================================================
# 7. INTERVAL SCORE
# =============================================================================

interval_score <- function(
    lower,
    upper,
    y,
    alpha
) {
  
  (
    upper -
      lower
  ) +
    
    (
      2 /
        alpha
    ) *
    (
      lower -
        y
    ) *
    as.numeric(
      y <
        lower
    ) +
    
    (
      2 /
        alpha
    ) *
    (
      y -
        upper
    ) *
    as.numeric(
      y >
        upper
    )
}



# =============================================================================
# 8. STANDARD WIS
# =============================================================================

calculate_wis <- function(
    y,
    q05,
    q10,
    q25,
    q50,
    q75,
    q90,
    q95
) {
  
  absolute_error <-
    abs(
      y -
        q50
    )
  
  
  IS50 <-
    interval_score(
      q25,
      q75,
      y,
      0.50
    )
  
  
  IS80 <-
    interval_score(
      q10,
      q90,
      y,
      0.20
    )
  
  
  IS90 <-
    interval_score(
      q05,
      q95,
      y,
      0.10
    )
  
  
  (
    0.5 *
      absolute_error +
      
      0.25 *
      IS50 +
      
      0.10 *
      IS80 +
      
      0.05 *
      IS90
  ) /
    3.5
}



# =============================================================================
# 9. EMPIRICAL CRPS
# =============================================================================

empirical_crps <- function(
    draws,
    y
) {
  
  draws <-
    draws[
      is.finite(
        draws
      )
    ]
  
  
  draws <-
    sort(
      draws
    )
  
  
  D <-
    length(
      draws
    )
  
  
  if (
    D <
    2L
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  term1 <-
    mean(
      abs(
        draws -
          y
      )
    )
  
  
  idx <-
    seq_len(
      D
    )
  
  
  term2 <-
    sum(
      (
        2 *
          idx -
          D -
          1
      ) *
        draws
    ) /
    D^2
  
  
  term1 -
    term2
}



# =============================================================================
# 10. READ DATA
# =============================================================================

if (
  !file.exists(
    STATE_DATA_FILE
  )
) {
  
  stop(
    paste0(
      "Cannot find data file: ",
      STATE_DATA_FILE
    )
  )
}


state_raw <-
  readr::read_csv(
    STATE_DATA_FILE,
    show_col_types = FALSE
  )



# =============================================================================
# 11. REQUIRED COLUMNS
# =============================================================================

required_columns <- c(
  "season",
  "state",
  "population",
  "epi_week_of_season",
  "cases"
)


missing_columns <-
  setdiff(
    required_columns,
    names(
      state_raw
    )
  )


if (
  length(
    missing_columns
  ) >
  0L
) {
  
  stop(
    paste0(
      "Missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}



# =============================================================================
# 12. BASIC DATA PREPARATION
# =============================================================================

dat <-
  state_raw %>%
  
  dplyr::mutate(
    
    season =
      normalize_season(
        .data$season
      ),
    
    season_index =
      match(
        .data$season,
        SEASON_LEVELS
      ),
    
    week =
      as.integer(
        .data$epi_week_of_season
      ),
    
    population =
      as.numeric(
        .data$population
      ),
    
    cases =
      as.numeric(
        .data$cases
      )
  ) %>%
  
  dplyr::filter(
    
    !is.na(
      .data$season
    ),
    
    !is.na(
      .data$state
    ),
    
    !is.na(
      .data$week
    ),
    
    !is.na(
      .data$population
    ),
    
    .data$population >
      0,
    
    !is.na(
      .data$cases
    )
  ) %>%
  
  dplyr::arrange(
    .data$season_index,
    .data$state,
    .data$week
  )



# =============================================================================
# 13. STATE INDEX
# =============================================================================

STATE_LEVELS <-
  sort(
    unique(
      dat$state
    )
  )


N_STATES <-
  length(
    STATE_LEVELS
  )


dat <-
  dat %>%
  
  dplyr::mutate(
    
    state_id =
      match(
        .data$state,
        STATE_LEVELS
      )
  )


cat(
  "\nNumber of states/FCT: ",
  N_STATES,
  "\n",
  sep = ""
)


if (
  N_STATES !=
  37L
) {
  
  warning(
    paste0(
      "Expected 37 states/FCT; found ",
      N_STATES
    )
  )
}



# =============================================================================
# 14. INCIDENCE
# =============================================================================

dat <-
  dat %>%
  
  dplyr::mutate(
    
    incidence =
      100000 *
      .data$cases /
      .data$population,
    
    log_incidence =
      log1p(
        .data$incidence
      )
  )



# =============================================================================
# 15. NATIONAL POPULATION-WEIGHTED EPIDEMIC LEVEL
# =============================================================================

national_data <-
  dat %>%
  
  dplyr::group_by(
    .data$season,
    .data$season_index,
    .data$week
  ) %>%
  
  dplyr::summarise(
    
    national_cases =
      sum(
        .data$cases,
        na.rm = TRUE
      ),
    
    national_population =
      sum(
        .data$population,
        na.rm = TRUE
      ),
    
    national_incidence =
      100000 *
      .data$national_cases /
      .data$national_population,
    
    national_log_incidence =
      log1p(
        .data$national_incidence
      ),
    
    .groups =
      "drop"
  )


dat <-
  dat %>%
  
  dplyr::left_join(
    
    national_data,
    
    by =
      c(
        "season",
        "season_index",
        "week"
      )
  )



# =============================================================================
# 16. MEMORY
#
# memory_lag1 = current X_t
# memory_lag2 = X_(t-1)
# ...
# memory_lag8 = X_(t-7)
# =============================================================================

dat <-
  dat %>%
  
  dplyr::group_by(
    .data$season,
    .data$state
  ) %>%
  
  dplyr::arrange(
    .data$week,
    .by_group = TRUE
  )


for (
  lag_index in
  0:(L_MEMORY - 1L)
) {
  
  variable_name <-
    paste0(
      "memory_lag",
      lag_index + 1L
    )
  
  
  dat <-
    dat %>%
    
    dplyr::mutate(
      
      !!variable_name :=
        dplyr::lag(
          .data$log_incidence,
          lag_index
        )
    )
}


dat <-
  dat %>%
  dplyr::ungroup()


memory_names <-
  paste0(
    "memory_lag",
    seq_len(
      L_MEMORY
    )
  )


memory_matrix <-
  as.matrix(
    dat[
      ,
      memory_names
    ]
  )


dat$memory <-
  as.numeric(
    memory_matrix %*%
      MEMORY_WEIGHTS
  )



# =============================================================================
# 17. NEXT-WEEK OUTCOME AND CAUSAL GROWTH
# =============================================================================

dat <-
  dat %>%
  
  dplyr::group_by(
    .data$season,
    .data$state
  ) %>%
  
  dplyr::arrange(
    .data$week,
    .by_group = TRUE
  ) %>%
  
  dplyr::mutate(
    
    cases_next =
      dplyr::lead(
        .data$cases,
        1
      ),
    
    log_incidence_next =
      dplyr::lead(
        .data$log_incidence,
        1
      ),
    
    growth_raw =
      .data$log_incidence -
      dplyr::lag(
        .data$log_incidence,
        1
      )
  ) %>%
  
  dplyr::ungroup()



# =============================================================================
# 18. TRAINING DATA
# =============================================================================

get_training_rows <- function(
    target_season_index,
    origin
) {
  
  dat %>%
    
    dplyr::filter(
      
      !is.na(
        .data$cases_next
      ),
      
      !is.na(
        .data$log_incidence_next
      ),
      
      !is.na(
        .data$memory
      ),
      
      !is.na(
        .data$growth_raw
      ),
      
      (
        .data$season_index <
          target_season_index
      ) |
        
        (
          .data$season_index ==
            target_season_index &
            
            .data$week <=
            origin -
            1L
        )
    )
}



# =============================================================================
# 19. FORECAST ROW
# =============================================================================

get_forecast_rows <- function(
    target_season_index,
    origin
) {
  
  dat %>%
    
    dplyr::filter(
      
      .data$season_index ==
        target_season_index,
      
      .data$week ==
        origin,
      
      !is.na(
        .data$cases_next
      ),
      
      !is.na(
        .data$memory
      ),
      
      !is.na(
        .data$growth_raw
      )
    ) %>%
    
    dplyr::arrange(
      .data$state_id
    )
}



# =============================================================================
# 20. CONNECTIVITY ESTIMATION
# =============================================================================

estimate_connectivity <- function(
    train_data,
    growth_scale
) {
  
  train_local <-
    train_data %>%
    
    dplyr::mutate(
      
      p_growth =
        plogis(
          .data$growth_raw /
            growth_scale
        )
    )
  
  
  edge_results <-
    list()
  
  
  counter <-
    1L
  
  
  for (
    target_idx in
    seq_len(
      N_STATES
    )
  ) {
    
    target_data <-
      train_local %>%
      
      dplyr::filter(
        .data$state_id ==
          target_idx
      ) %>%
      
      dplyr::transmute(
        
        season =
          .data$season,
        
        week =
          .data$week,
        
        y_next =
          .data$log_incidence_next,
        
        own_memory =
          .data$memory,
        
        national =
          .data$national_log_incidence,
        
        p_growth =
          .data$p_growth
      )
    
    
    for (
      source_idx in
      seq_len(
        N_STATES
      )
    ) {
      
      if (
        source_idx ==
        target_idx
      ) {
        
        next
      }
      
      
      source_data <-
        train_local %>%
        
        dplyr::filter(
          .data$state_id ==
            source_idx
        ) %>%
        
        dplyr::transmute(
          
          season =
            .data$season,
          
          week =
            .data$week,
          
          source_current =
            .data$log_incidence
        )
      
      
      edge_data <-
        target_data %>%
        
        dplyr::inner_join(
          
          source_data,
          
          by =
            c(
              "season",
              "week"
            )
        ) %>%
        
        dplyr::filter(
          
          stats::complete.cases(
            
            .data$y_next,
            
            .data$own_memory,
            
            .data$national,
            
            .data$p_growth,
            
            .data$source_current
          )
        )
      
      
      n_edge <-
        nrow(
          edge_data
        )
      
      
      source_coefficient <-
        NA_real_
      
      
      if (
        n_edge >=
        MIN_EDGE_N
      ) {
        
        edge_fit <-
          try(
            
            stats::lm(
              
              y_next ~
                own_memory +
                national +
                p_growth +
                source_current,
              
              data =
                edge_data
            ),
            
            silent =
              TRUE
          )
        
        
        if (
          !inherits(
            edge_fit,
            "try-error"
          )
        ) {
          
          coefficient_vector <-
            stats::coef(
              edge_fit
            )
          
          
          if (
            "source_current" %in%
            names(
              coefficient_vector
            )
          ) {
            
            source_coefficient <-
              unname(
                coefficient_vector[
                  "source_current"
                ]
              )
          }
        }
      }
      
      
      edge_results[[
        counter
      ]] <-
        
        tibble::tibble(
          
          target_id =
            target_idx,
          
          source_id =
            source_idx,
          
          n =
            n_edge,
          
          coefficient =
            source_coefficient
        )
      
      
      counter <-
        counter +
        1L
    }
  }
  
  
  edges <-
    dplyr::bind_rows(
      edge_results
    )
  
  
  estimable_edges <-
    sum(
      is.finite(
        edges$coefficient
      )
    )
  
  
  positive_edges <-
    sum(
      
      is.finite(
        edges$coefficient
      ) &
        
        edges$coefficient >
        0
    )
  
  
  W <-
    matrix(
      0,
      nrow =
        N_STATES,
      ncol =
        N_STATES
    )
  
  
  for (
    target_idx in
    seq_len(
      N_STATES
    )
  ) {
    
    target_edges <-
      edges %>%
      
      dplyr::filter(
        .data$target_id ==
          target_idx
      )
    
    
    weights <-
      ifelse(
        
        is.finite(
          target_edges$coefficient
        ),
        
        pmax(
          target_edges$coefficient,
          0
        ),
        
        0
      )
    
    
    if (
      sum(
        weights
      ) >
      0
    ) {
      
      weights <-
        weights /
        sum(
          weights
        )
      
    } else {
      
      weights[] <-
        1 /
        (
          N_STATES -
            1L
        )
    }
    
    
    W[
      target_idx,
      target_edges$source_id
    ] <-
      weights
  }
  
  
  diag(W) <-
    0
  
  
  row_error <-
    max(
      abs(
        rowSums(W) -
          1
      )
    )
  
  
  if (
    row_error >
    1e-8
  ) {
    
    stop(
      "Connectivity rows do not sum to one."
    )
  }
  
  
  list(
    
    W =
      W,
    
    edges =
      edges,
    
    estimable_edges =
      estimable_edges,
    
    positive_edges =
      positive_edges
  )
}



# =============================================================================
# 21. BORROWING
# =============================================================================

calculate_borrowing <- function(
    data_rows,
    W
) {
  
  result <-
    rep(
      NA_real_,
      nrow(
        data_rows
      )
    )
  
  
  blocks <-
    data_rows %>%
    
    dplyr::distinct(
      .data$season,
      .data$week
    )
  
  
  for (
    block_idx in
    seq_len(
      nrow(
        blocks
      )
    )
  ) {
    
    block_season <-
      blocks$season[
        block_idx
      ]
    
    
    block_week <-
      blocks$week[
        block_idx
      ]
    
    
    indices <-
      which(
        
        data_rows$season ==
          block_season &
          
          data_rows$week ==
          block_week
      )
    
    
    block <-
      data_rows[
        indices,
        ,
        drop = FALSE
      ] %>%
      
      dplyr::arrange(
        .data$state_id
      )
    
    
    if (
      nrow(
        block
      ) !=
      N_STATES
    ) {
      
      next
    }
    
    
    if (
      !all(
        block$state_id ==
        seq_len(
          N_STATES
        )
      )
    ) {
      
      stop(
        "State ordering error."
      )
    }
    
    
    source_signal <-
      block$log_incidence
    
    
    borrowed_signal <-
      as.numeric(
        W %*%
          source_signal
      )
    
    
    original_state_ids <-
      data_rows$state_id[
        indices
      ]
    
    
    result[
      indices
    ] <-
      
      borrowed_signal[
        original_state_ids
      ]
  }
  
  
  result
}



# =============================================================================
# 22. PREPARE ONE ORIGIN
# =============================================================================

prepare_origin <- function(
    target_season_index,
    origin
) {
  
  target_season <-
    SEASON_LEVELS[
      target_season_index
    ]
  
  
  cat(
    "\nPREPARING ",
    target_season,
    " ORIGIN ",
    origin,
    "\n",
    sep = ""
  )
  
  
  train <-
    get_training_rows(
      target_season_index,
      origin
    )
  
  
  pred <-
    get_forecast_rows(
      target_season_index,
      origin
    )
  
  
  if (
    nrow(
      pred
    ) !=
    N_STATES
  ) {
    
    stop(
      paste0(
        "Expected ",
        N_STATES,
        " prediction rows; obtained ",
        nrow(
          pred
        )
      )
    )
  }
  
  
  # ===========================================================================
  # TRAIN-ONLY PHASE SCALE
  # ===========================================================================
  
  growth_scale <-
    safe_sd(
      train$growth_raw
    )
  
  
  train <-
    train %>%
    
    dplyr::mutate(
      
      p_growth =
        plogis(
          .data$growth_raw /
            growth_scale
        )
    )
  
  
  pred <-
    pred %>%
    
    dplyr::mutate(
      
      p_growth =
        plogis(
          .data$growth_raw /
            growth_scale
        )
    )
  
  
  # ===========================================================================
  # CONNECTIVITY
  # ===========================================================================
  
  connectivity <-
    estimate_connectivity(
      train,
      growth_scale
    )
  
  
  W <-
    connectivity$W
  
  
  train$borrowing <-
    calculate_borrowing(
      train,
      W
    )
  
  
  pred$borrowing <-
    calculate_borrowing(
      pred,
      W
    )
  
  
  # ===========================================================================
  # FINAL FEATURES
  # ===========================================================================
  
  train <-
    train %>%
    
    dplyr::mutate(
      
      borrow_phase =
        .data$borrowing *
        .data$p_growth,
      
      log_offset =
        log(
          .data$population /
            100000
        )
    )
  
  
  pred <-
    pred %>%
    
    dplyr::mutate(
      
      borrow_phase =
        .data$borrowing *
        .data$p_growth,
      
      log_offset =
        log(
          .data$population /
            100000
        )
    )
  
  
  # ===========================================================================
  # COMPLETE CASES
  # ===========================================================================
  
  train <-
    train %>%
    
    dplyr::filter(
      
      stats::complete.cases(
        
        .data$cases_next,
        
        .data$population,
        
        .data$memory,
        
        .data$national_log_incidence,
        
        .data$borrowing,
        
        .data$p_growth,
        
        .data$borrow_phase,
        
        .data$log_offset
      )
    )
  
  
  pred <-
    pred %>%
    
    dplyr::filter(
      
      stats::complete.cases(
        
        .data$cases_next,
        
        .data$population,
        
        .data$memory,
        
        .data$national_log_incidence,
        
        .data$borrowing,
        
        .data$p_growth,
        
        .data$borrow_phase,
        
        .data$log_offset
      )
    )
  
  
  if (
    nrow(
      pred
    ) !=
    N_STATES
  ) {
    
    stop(
      "Prediction rows incomplete after feature generation."
    )
  }
  
  
  # ===========================================================================
  # CAUSAL ORIGIN INTENSITY FOR DIAGNOSTICS ONLY
  #
  # This is NOT used in M4.
  # ===========================================================================
  
  intensity_center <-
    mean(
      train$log_incidence,
      na.rm = TRUE
    )
  
  
  intensity_scale <-
    safe_sd(
      train$log_incidence
    )
  
  
  pred <-
    pred %>%
    
    dplyr::mutate(
      
      origin_intensity_z =
        (
          .data$log_incidence -
            intensity_center
        ) /
        intensity_scale
    )
  
  
  # ===========================================================================
  # FINITE CHECK
  # ===========================================================================
  
  check_variables <- c(
    
    "memory",
    
    "national_log_incidence",
    
    "borrowing",
    
    "p_growth",
    
    "borrow_phase",
    
    "log_offset"
  )
  
  
  for (
    variable_name in
    check_variables
  ) {
    
    check_finite(
      
      train[[
        variable_name
      ]],
      
      paste0(
        "train ",
        variable_name
      )
    )
    
    
    check_finite(
      
      pred[[
        variable_name
      ]],
      
      paste0(
        "prediction ",
        variable_name
      )
    )
  }
  
  
  list(
    
    season =
      target_season,
    
    season_index =
      target_season_index,
    
    origin =
      origin,
    
    target_week =
      origin +
      1L,
    
    train =
      train,
    
    pred =
      pred,
    
    growth_scale =
      growth_scale,
    
    intensity_center =
      intensity_center,
    
    intensity_scale =
      intensity_scale,
    
    connectivity =
      connectivity
  )
}



# =============================================================================
# 23. STAN DATA
# =============================================================================

build_stan_data <- function(
    prepared
) {
  
  train <-
    prepared$train
  
  
  pred <-
    prepared$pred
  
  
  list(
    
    N_train =
      nrow(
        train
      ),
    
    N_pred =
      nrow(
        pred
      ),
    
    J =
      N_STATES,
    
    
    y =
      as.integer(
        train$cases_next
      ),
    
    
    state_id =
      as.integer(
        train$state_id
      ),
    
    
    state_id_pred =
      as.integer(
        pred$state_id
      ),
    
    
    log_offset =
      as.vector(
        train$log_offset
      ),
    
    
    log_offset_pred =
      as.vector(
        pred$log_offset
      ),
    
    
    memory =
      as.vector(
        train$memory
      ),
    
    
    memory_pred =
      as.vector(
        pred$memory
      ),
    
    
    national =
      as.vector(
        train$national_log_incidence
      ),
    
    
    national_pred =
      as.vector(
        pred$national_log_incidence
      ),
    
    
    borrowing =
      as.vector(
        train$borrowing
      ),
    
    
    borrowing_pred =
      as.vector(
        pred$borrowing
      ),
    
    
    phase =
      as.vector(
        train$p_growth
      ),
    
    
    phase_pred =
      as.vector(
        pred$p_growth
      ),
    
    
    borrow_phase =
      as.vector(
        train$borrow_phase
      ),
    
    
    borrow_phase_pred =
      as.vector(
        pred$borrow_phase
      )
  )
}



# =============================================================================
# 24. FINAL M4 STAN MODEL
# =============================================================================

stan_code_m4 <- '

data {

  int<lower=1> N_train;

  int<lower=1> N_pred;

  int<lower=1> J;


  array[N_train]
    int<lower=0>
    y;


  array[N_train]
    int<lower=1,upper=J>
    state_id;


  array[N_pred]
    int<lower=1,upper=J>
    state_id_pred;


  vector[N_train]
    log_offset;

  vector[N_pred]
    log_offset_pred;


  vector[N_train]
    memory;

  vector[N_pred]
    memory_pred;


  vector[N_train]
    national;

  vector[N_pred]
    national_pred;


  vector[N_train]
    borrowing;

  vector[N_pred]
    borrowing_pred;


  vector[N_train]
    phase;

  vector[N_pred]
    phase_pred;


  vector[N_train]
    borrow_phase;

  vector[N_pred]
    borrow_phase_pred;
}


parameters {

  real alpha;


  vector[J]
    state_raw;


  real<lower=0>
    sigma_state;


  real beta_memory;

  real beta_national;

  real beta_borrow;

  real beta_phase;

  real beta_interaction;


  real log_phi;
}


transformed parameters {

  vector[J]
    alpha_state;


  alpha_state =
    sigma_state *
    state_raw;
}


model {

  vector[N_train]
    eta;


  real phi;


  alpha ~
    normal(
      0,
      2
    );


  state_raw ~
    normal(
      0,
      1
    );


  sigma_state ~
    normal(
      0,
      1
    );


  beta_memory ~
    normal(
      0,
      1
    );


  beta_national ~
    normal(
      0,
      1
    );


  beta_borrow ~
    normal(
      0,
      1
    );


  beta_phase ~
    normal(
      0,
      1
    );


  beta_interaction ~
    normal(
      0,
      1
    );


  log_phi ~
    normal(
      log(5),
      0.5
    );


  phi =
    0.05 +
    exp(
      log_phi
    );


  eta =

      log_offset

      + alpha

      + alpha_state[
          state_id
        ]

      + beta_memory *
        memory

      + beta_national *
        national

      + beta_borrow *
        borrowing

      + beta_phase *
        phase

      + beta_interaction *
        borrow_phase;


  y ~
    neg_binomial_2_log(
      eta,
      phi
    );
}


generated quantities {

  vector[N_pred]
    eta_pred;


  real phi_pred;


  array[N_pred]
    int y_rep;


  eta_pred =

      log_offset_pred

      + alpha

      + alpha_state[
          state_id_pred
        ]

      + beta_memory *
        memory_pred

      + beta_national *
        national_pred

      + beta_borrow *
        borrowing_pred

      + beta_phase *
        phase_pred

      + beta_interaction *
        borrow_phase_pred;


  phi_pred =
    0.05 +
    exp(
      log_phi
    );


  for (
    n in
    1:N_pred
  ) {

    y_rep[n] =
      neg_binomial_2_log_rng(
        eta_pred[n],
        phi_pred
      );
  }
}
'



# =============================================================================
# 25. WRITE AND COMPILE STAN MODEL
# =============================================================================

M4_STAN_FILE <-
  file.path(
    STAN_DIR,
    "PACE_FLU_FINAL_M4.stan"
  )


writeLines(
  stan_code_m4,
  M4_STAN_FILE
)


cat(
  "\nCompiling final M4...\n"
)


model_m4 <-
  cmdstanr::cmdstan_model(
    M4_STAN_FILE
  )



# =============================================================================
# 26. PARAMETER NAMES
# =============================================================================

PARAMETERS <- c(
  
  "alpha",
  
  "sigma_state",
  
  "beta_memory",
  
  "beta_national",
  
  "beta_borrow",
  
  "beta_phase",
  
  "beta_interaction",
  
  "log_phi"
)


CONVERGENCE_PARAMETERS <-
  c(
    "state_raw",
    PARAMETERS
  )



# =============================================================================
# 27. PARAMETER SUMMARY
# =============================================================================

extract_parameter_summary <- function(
    fit,
    season,
    origin
) {
  
  draws_matrix <-
    fit$draws(
      variables =
        PARAMETERS,
      format =
        "matrix"
    )
  
  
  dplyr::bind_rows(
    
    lapply(
      
      PARAMETERS,
      
      function(
    variable_name
      ) {
        
        x <-
          as.numeric(
            draws_matrix[
              ,
              variable_name
            ]
          )
        
        
        tibble::tibble(
          
          season =
            season,
          
          origin =
            origin,
          
          variable =
            variable_name,
          
          mean =
            mean(
              x,
              na.rm = TRUE
            ),
          
          median =
            stats::median(
              x,
              na.rm = TRUE
            ),
          
          sd =
            stats::sd(
              x,
              na.rm = TRUE
            ),
          
          q05 =
            as.numeric(
              stats::quantile(
                x,
                0.05,
                na.rm = TRUE,
                names = FALSE
              )
            ),
          
          q25 =
            as.numeric(
              stats::quantile(
                x,
                0.25,
                na.rm = TRUE,
                names = FALSE
              )
            ),
          
          q75 =
            as.numeric(
              stats::quantile(
                x,
                0.75,
                na.rm = TRUE,
                names = FALSE
              )
            ),
          
          q95 =
            as.numeric(
              stats::quantile(
                x,
                0.95,
                na.rm = TRUE,
                names = FALSE
              )
            ),
          
          probability_positive =
            mean(
              x >
                0,
              na.rm = TRUE
            )
        )
      }
    )
  )
}



# =============================================================================
# 28. CONVERGENCE SUMMARY
# =============================================================================

extract_convergence <- function(
    fit,
    season,
    origin
) {
  
  sm <-
    fit$summary(
      variables =
        CONVERGENCE_PARAMETERS
    )
  
  
  finite_rhat <-
    sm$rhat[
      is.finite(
        sm$rhat
      )
    ]
  
  
  finite_bulk <-
    sm$ess_bulk[
      is.finite(
        sm$ess_bulk
      )
    ]
  
  
  finite_tail <-
    sm$ess_tail[
      is.finite(
        sm$ess_tail
      )
    ]
  
  
  max_rhat <-
    if (
      length(
        finite_rhat
      ) >
      0L
    ) {
      
      max(
        finite_rhat
      )
      
    } else {
      
      NA_real_
    }
  
  
  min_bulk <-
    if (
      length(
        finite_bulk
      ) >
      0L
    ) {
      
      min(
        finite_bulk
      )
      
    } else {
      
      NA_real_
    }
  
  
  min_tail <-
    if (
      length(
        finite_tail
      ) >
      0L
    ) {
      
      min(
        finite_tail
      )
      
    } else {
      
      NA_real_
    }
  
  
  diagnostic <-
    fit$diagnostic_summary()
  
  
  divergences <-
    sum(
      diagnostic$num_divergent,
      na.rm = TRUE
    )
  
  
  treedepth_hits <-
    sum(
      diagnostic$num_max_treedepth,
      na.rm = TRUE
    )
  
  
  convergence_ok <-
    
    is.finite(
      max_rhat
    ) &&
    
    max_rhat <
    1.01 &&
    
    is.finite(
      min_bulk
    ) &&
    
    min_bulk >
    400 &&
    
    is.finite(
      min_tail
    ) &&
    
    min_tail >
    400 &&
    
    divergences ==
    0 &&
    
    treedepth_hits ==
    0
  
  
  tibble::tibble(
    
    season =
      season,
    
    origin =
      origin,
    
    max_rhat =
      max_rhat,
    
    min_ess_bulk =
      min_bulk,
    
    min_ess_tail =
      min_tail,
    
    divergences =
      divergences,
    
    treedepth_hits =
      treedepth_hits,
    
    convergence_ok =
      convergence_ok
  )
}



# =============================================================================
# 29. SCORE FORECAST
# =============================================================================

score_fit <- function(
    fit,
    prepared
) {
  
  pred <-
    prepared$pred
  
  
  draws <-
    fit$draws(
      variables =
        "y_rep",
      format =
        "matrix"
    )
  
  
  expected_names <-
    paste0(
      "y_rep[",
      seq_len(
        nrow(
          pred
        )
      ),
      "]"
    )
  
  
  draws <-
    draws[
      ,
      expected_names,
      drop = FALSE
    ]
  
  
  if (
    nrow(
      draws
    ) >
    MAX_PRED_DRAWS
  ) {
    
    keep <-
      unique(
        round(
          seq(
            1,
            nrow(
              draws
            ),
            length.out =
              MAX_PRED_DRAWS
          )
        )
      )
    
    
    draws <-
      draws[
        keep,
        ,
        drop = FALSE
      ]
  }
  
  
  result <-
    vector(
      "list",
      nrow(
        pred
      )
    )
  
  
  probs <- c(
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95
  )
  
  
  for (
    row_idx in
    seq_len(
      nrow(
        pred
      )
    )
  ) {
    
    draws_i <-
      as.numeric(
        draws[
          ,
          row_idx
        ]
      )
    
    
    q <-
      stats::quantile(
        
        draws_i,
        
        probs =
          probs,
        
        na.rm =
          TRUE,
        
        names =
          FALSE,
        
        type =
          8
      )
    
    
    actual <-
      pred$cases_next[
        row_idx
      ]
    
    
    result[[
      row_idx
    ]] <-
      
      tibble::tibble(
        
        season =
          prepared$season,
        
        origin =
          prepared$origin,
        
        target_week =
          prepared$target_week,
        
        state =
          pred$state[
            row_idx
          ],
        
        state_id =
          pred$state_id[
            row_idx
          ],
        
        population =
          pred$population[
            row_idx
          ],
        
        cases_at_origin =
          pred$cases[
            row_idx
          ],
        
        incidence_at_origin =
          pred$incidence[
            row_idx
          ],
        
        log_incidence_at_origin =
          pred$log_incidence[
            row_idx
          ],
        
        origin_intensity_z =
          pred$origin_intensity_z[
            row_idx
          ],
        
        p_growth =
          pred$p_growth[
            row_idx
          ],
        
        borrowing =
          pred$borrowing[
            row_idx
          ],
        
        borrow_phase =
          pred$borrow_phase[
            row_idx
          ],
        
        actual =
          actual,
        
        q05 =
          q[1],
        
        q10 =
          q[2],
        
        q25 =
          q[3],
        
        q50 =
          q[4],
        
        q75 =
          q[5],
        
        q90 =
          q[6],
        
        q95 =
          q[7],
        
        error =
          q[4] -
          actual,
        
        underprediction_error =
          actual -
          q[4],
        
        absolute_error =
          abs(
            q[4] -
              actual
          ),
        
        squared_error =
          (
            q[4] -
              actual
          )^2,
        
        WIS =
          calculate_wis(
            
            actual,
            
            q[1],
            
            q[2],
            
            q[3],
            
            q[4],
            
            q[5],
            
            q[6],
            
            q[7]
          ),
        
        CRPS =
          empirical_crps(
            draws_i,
            actual
          ),
        
        covered50 =
          actual >=
          q[3] &&
          actual <=
          q[5],
        
        covered80 =
          actual >=
          q[2] &&
          actual <=
          q[6],
        
        covered90 =
          actual >=
          q[1] &&
          actual <=
          q[7],
        
        above90 =
          actual >
          q[7],
        
        below90 =
          actual <
          q[1],
        
        width50 =
          q[5] -
          q[3],
        
        width80 =
          q[6] -
          q[2],
        
        width90 =
          q[7] -
          q[1]
      )
  }
  
  
  dplyr::bind_rows(
    result
  )
}



# =============================================================================
# 30. STORAGE
# =============================================================================

feature_results <-
  list()


parameter_results <-
  list()


convergence_results <-
  list()


forecast_results <-
  list()



# =============================================================================
# 31. FINAL 24-ORIGIN VALIDATION
# =============================================================================

for (
  grid_idx in
  seq_len(
    nrow(
      VALIDATION_GRID
    )
  )
) {
  
  TARGET_SEASON <-
    VALIDATION_GRID$season[
      grid_idx
    ]
  
  
  TARGET_SEASON_INDEX <-
    VALIDATION_GRID$season_index[
      grid_idx
    ]
  
  
  TARGET_ORIGIN <-
    VALIDATION_GRID$origin[
      grid_idx
    ]
  
  
  TARGET_WEEK <-
    TARGET_ORIGIN +
    1L
  
  
  cat(
    
    "\n\n============================================================\n",
    
    "FINAL M4: ",
    
    TARGET_SEASON,
    
    " origin ",
    
    TARGET_ORIGIN,
    
    " -> ",
    
    TARGET_WEEK,
    
    "\n============================================================\n",
    
    sep = ""
  )
  
  
  # ===========================================================================
  # PREPARE
  # ===========================================================================
  
  prepared <-
    prepare_origin(
      
      TARGET_SEASON_INDEX,
      
      TARGET_ORIGIN
    )
  
  
  # ===========================================================================
  # FEATURE SUMMARY
  # ===========================================================================
  
  feature_results[[
    grid_idx
  ]] <-
    
    tibble::tibble(
      
      season =
        TARGET_SEASON,
      
      origin =
        TARGET_ORIGIN,
      
      target_week =
        TARGET_WEEK,
      
      training_rows =
        nrow(
          prepared$train
        ),
      
      prediction_rows =
        nrow(
          prepared$pred
        ),
      
      estimable_edges =
        prepared$connectivity$estimable_edges,
      
      positive_edges =
        prepared$connectivity$positive_edges,
      
      learned_network =
        prepared$connectivity$estimable_edges >
        0,
      
      growth_scale =
        prepared$growth_scale,
      
      intensity_center =
        prepared$intensity_center,
      
      intensity_scale =
        prepared$intensity_scale
    )
  
  
  # ===========================================================================
  # STAN DATA
  # ===========================================================================
  
  stan_data <-
    build_stan_data(
      prepared
    )
  
  
  # ===========================================================================
  # FIT
  # ===========================================================================
  
  cat(
    "\nFITTING FINAL M4\n"
  )
  
  
  fit <-
    model_m4$sample(
      
      data =
        stan_data,
      
      seed =
        SEED +
        grid_idx *
        100L,
      
      chains =
        N_CHAINS,
      
      parallel_chains =
        N_PARALLEL_CHAINS,
      
      iter_warmup =
        N_WARMUP,
      
      iter_sampling =
        N_SAMPLING,
      
      init =
        STAN_INIT,
      
      step_size =
        INITIAL_STEP_SIZE,
      
      adapt_delta =
        ADAPT_DELTA,
      
      max_treedepth =
        MAX_TREEDEPTH,
      
      refresh =
        100
    )
  
  
  # ===========================================================================
  # SAVE FIT IMMEDIATELY
  # ===========================================================================
  
  origin_tag <-
    paste0(
      
      gsub(
        "/",
        "_",
        TARGET_SEASON
      ),
      
      "_origin_",
      
      TARGET_ORIGIN
    )
  
  
  saveRDS(
    
    fit,
    
    file.path(
      
      FIT_DIR,
      
      paste0(
        "M4_",
        origin_tag,
        ".rds"
      )
    )
  )
  
  
  # ===========================================================================
  # PARAMETERS
  # ===========================================================================
  
  parameter_results[[
    grid_idx
  ]] <-
    
    extract_parameter_summary(
      
      fit,
      
      TARGET_SEASON,
      
      TARGET_ORIGIN
    )
  
  
  # ===========================================================================
  # CONVERGENCE
  # ===========================================================================
  
  convergence_results[[
    grid_idx
  ]] <-
    
    extract_convergence(
      
      fit,
      
      TARGET_SEASON,
      
      TARGET_ORIGIN
    )
  
  
  # ===========================================================================
  # FORECASTS
  # ===========================================================================
  
  forecast_results[[
    grid_idx
  ]] <-
    
    score_fit(
      fit,
      prepared
    )
  
  
  # ===========================================================================
  # CHECKPOINT
  # ===========================================================================
  
  readr::write_csv(
    
    feature_results[[
      grid_idx
    ]],
    
    file.path(
      
      CHECKPOINT_DIR,
      
      paste0(
        "feature_",
        origin_tag,
        ".csv"
      )
    )
  )
  
  
  readr::write_csv(
    
    parameter_results[[
      grid_idx
    ]],
    
    file.path(
      
      CHECKPOINT_DIR,
      
      paste0(
        "parameters_",
        origin_tag,
        ".csv"
      )
    )
  )
  
  
  readr::write_csv(
    
    convergence_results[[
      grid_idx
    ]],
    
    file.path(
      
      CHECKPOINT_DIR,
      
      paste0(
        "convergence_",
        origin_tag,
        ".csv"
      )
    )
  )
  
  
  readr::write_csv(
    
    forecast_results[[
      grid_idx
    ]],
    
    file.path(
      
      CHECKPOINT_DIR,
      
      paste0(
        "forecast_",
        origin_tag,
        ".csv"
      )
    )
  )
  
  
  cat(
    "\nCONVERGENCE\n"
  )
  
  
  safe_print(
    convergence_results[[
      grid_idx
    ]]
  )
  
  
  rm(
    fit
  )
  
  
  invisible(
    gc()
  )
}



# =============================================================================
# 32. COMBINE FINAL RESULTS
# =============================================================================

feature_summary <-
  dplyr::bind_rows(
    feature_results
  )


parameter_summary <-
  dplyr::bind_rows(
    parameter_results
  )


convergence_summary <-
  dplyr::bind_rows(
    convergence_results
  )


all_forecasts <-
  dplyr::bind_rows(
    forecast_results
  )



# =============================================================================
# 33. PAIRING / COMPLETENESS CHECK
# =============================================================================

pairing_check <-
  VALIDATION_GRID %>%
  
  dplyr::left_join(
    
    all_forecasts %>%
      
      dplyr::count(
        .data$season,
        .data$origin,
        name =
          "n_state_forecasts"
      ),
    
    by =
      c(
        "season",
        "origin"
      )
  ) %>%
  
  dplyr::mutate(
    
    expected_states =
      N_STATES,
    
    complete_origin =
      .data$n_state_forecasts ==
      N_STATES
  )



# =============================================================================
# 34. ORIGIN-LEVEL PERFORMANCE
# =============================================================================

origin_performance <-
  all_forecasts %>%
  
  dplyr::group_by(
    .data$season,
    .data$origin,
    .data$target_week
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error,
          na.rm = TRUE
        )
      ),
    
    bias =
      mean(
        .data$error,
        na.rm = TRUE
      ),
    
    median_bias =
      stats::median(
        .data$error,
        na.rm = TRUE
      ),
    
    mean_WIS =
      mean(
        .data$WIS,
        na.rm = TRUE
      ),
    
    median_WIS =
      stats::median(
        .data$WIS,
        na.rm = TRUE
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    median_CRPS =
      stats::median(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    coverage50 =
      mean(
        .data$covered50,
        na.rm = TRUE
      ),
    
    coverage80 =
      mean(
        .data$covered80,
        na.rm = TRUE
      ),
    
    coverage90 =
      mean(
        .data$covered90,
        na.rm = TRUE
      ),
    
    above90_rate =
      mean(
        .data$above90,
        na.rm = TRUE
      ),
    
    below90_rate =
      mean(
        .data$below90,
        na.rm = TRUE
      ),
    
    mean_width50 =
      mean(
        .data$width50,
        na.rm = TRUE
      ),
    
    mean_width80 =
      mean(
        .data$width80,
        na.rm = TRUE
      ),
    
    mean_width90 =
      mean(
        .data$width90,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )



# =============================================================================
# 35. OVERALL PERFORMANCE
#
# IMPORTANT:
# Median WIS and CRPS are calculated from the 888 state-origin scores.
# They are NOT copied from the means.
# =============================================================================

overall_performance <-
  all_forecasts %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    median_absolute_error =
      stats::median(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error,
          na.rm = TRUE
        )
      ),
    
    mean_bias =
      mean(
        .data$error,
        na.rm = TRUE
      ),
    
    median_bias =
      stats::median(
        .data$error,
        na.rm = TRUE
      ),
    
    mean_WIS =
      mean(
        .data$WIS,
        na.rm = TRUE
      ),
    
    median_WIS =
      stats::median(
        .data$WIS,
        na.rm = TRUE
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    median_CRPS =
      stats::median(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    coverage50 =
      mean(
        .data$covered50,
        na.rm = TRUE
      ),
    
    coverage80 =
      mean(
        .data$covered80,
        na.rm = TRUE
      ),
    
    coverage90 =
      mean(
        .data$covered90,
        na.rm = TRUE
      ),
    
    mean_width50 =
      mean(
        .data$width50,
        na.rm = TRUE
      ),
    
    mean_width80 =
      mean(
        .data$width80,
        na.rm = TRUE
      ),
    
    mean_width90 =
      mean(
        .data$width90,
        na.rm = TRUE
      )
  )



# =============================================================================
# 36. SEASON PERFORMANCE
# =============================================================================

season_performance <-
  all_forecasts %>%
  
  dplyr::group_by(
    .data$season
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error,
          na.rm = TRUE
        )
      ),
    
    mean_bias =
      mean(
        .data$error,
        na.rm = TRUE
      ),
    
    median_bias =
      stats::median(
        .data$error,
        na.rm = TRUE
      ),
    
    mean_WIS =
      mean(
        .data$WIS,
        na.rm = TRUE
      ),
    
    median_WIS =
      stats::median(
        .data$WIS,
        na.rm = TRUE
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    median_CRPS =
      stats::median(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    coverage50 =
      mean(
        .data$covered50,
        na.rm = TRUE
      ),
    
    coverage80 =
      mean(
        .data$covered80,
        na.rm = TRUE
      ),
    
    coverage90 =
      mean(
        .data$covered90,
        na.rm = TRUE
      ),
    
    mean_width90 =
      mean(
        .data$width90,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )



# =============================================================================
# 37. CALIBRATION SUMMARY
# =============================================================================

calibration_summary <-
  tibble::tibble(
    
    interval =
      c(
        "50%",
        "80%",
        "90%"
      ),
    
    nominal =
      c(
        0.50,
        0.80,
        0.90
      ),
    
    observed =
      c(
        
        mean(
          all_forecasts$covered50
        ),
        
        mean(
          all_forecasts$covered80
        ),
        
        mean(
          all_forecasts$covered90
        )
      )
  ) %>%
  
  dplyr::mutate(
    
    deviation =
      .data$observed -
      .data$nominal
  )



# =============================================================================
# 38. INTERVAL MISS DIRECTION
# =============================================================================

interval_miss_direction <-
  tibble::tibble(
    
    interval =
      c(
        "50%",
        "80%",
        "90%"
      ),
    
    above =
      c(
        
        sum(
          all_forecasts$actual >
            all_forecasts$q75
        ),
        
        sum(
          all_forecasts$actual >
            all_forecasts$q90
        ),
        
        sum(
          all_forecasts$actual >
            all_forecasts$q95
        )
      ),
    
    below =
      c(
        
        sum(
          all_forecasts$actual <
            all_forecasts$q25
        ),
        
        sum(
          all_forecasts$actual <
            all_forecasts$q10
        ),
        
        sum(
          all_forecasts$actual <
            all_forecasts$q05
        )
      ),
    
    covered =
      c(
        
        sum(
          all_forecasts$covered50
        ),
        
        sum(
          all_forecasts$covered80
        ),
        
        sum(
          all_forecasts$covered90
        )
      )
  ) %>%
  
  dplyr::mutate(
    
    total =
      .data$above +
      .data$below +
      .data$covered,
    
    proportion_above =
      .data$above /
      .data$total,
    
    proportion_below =
      .data$below /
      .data$total,
    
    proportion_covered =
      .data$covered /
      .data$total,
    
    proportion_misses_above =
      
      ifelse(
        
        (
          .data$above +
            .data$below
        ) >
          0,
        
        .data$above /
          (
            .data$above +
              .data$below
          ),
        
        NA_real_
      )
  )



# =============================================================================
# 39. CAUSAL ORIGIN-INTENSITY QUARTILES
#
# IMPORTANT:
# Quartiles are defined from information available at forecast origin.
# They are NOT based on future t+1 burden.
# =============================================================================

all_forecasts <-
  all_forecasts %>%
  
  dplyr::mutate(
    
    intensity_quartile =
      dplyr::ntile(
        .data$origin_intensity_z,
        4
      ),
    
    intensity_group =
      factor(
        
        .data$intensity_quartile,
        
        levels =
          1:4,
        
        labels =
          c(
            "Q1 Lowest",
            "Q2 Low-middle",
            "Q3 High-middle",
            "Q4 Highest"
          )
      )
  )



# =============================================================================
# 40. PERFORMANCE BY CAUSAL INTENSITY
# =============================================================================

intensity_performance <-
  all_forecasts %>%
  
  dplyr::group_by(
    .data$intensity_group
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    median_origin_cases =
      stats::median(
        .data$cases_at_origin,
        na.rm = TRUE
      ),
    
    mean_origin_cases =
      mean(
        .data$cases_at_origin,
        na.rm = TRUE
      ),
    
    median_origin_incidence =
      stats::median(
        .data$incidence_at_origin,
        na.rm = TRUE
      ),
    
    mean_actual =
      mean(
        .data$actual,
        na.rm = TRUE
      ),
    
    mean_forecast =
      mean(
        .data$q50,
        na.rm = TRUE
      ),
    
    mean_bias =
      mean(
        .data$error,
        na.rm = TRUE
      ),
    
    median_bias =
      stats::median(
        .data$error,
        na.rm = TRUE
      ),
    
    MAE =
      mean(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error,
          na.rm = TRUE
        )
      ),
    
    mean_WIS =
      mean(
        .data$WIS,
        na.rm = TRUE
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    coverage50 =
      mean(
        .data$covered50,
        na.rm = TRUE
      ),
    
    coverage80 =
      mean(
        .data$covered80,
        na.rm = TRUE
      ),
    
    coverage90 =
      mean(
        .data$covered90,
        na.rm = TRUE
      ),
    
    above90_rate =
      mean(
        .data$above90,
        na.rm = TRUE
      ),
    
    below90_rate =
      mean(
        .data$below90,
        na.rm = TRUE
      ),
    
    mean_width90 =
      mean(
        .data$width90,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )



# =============================================================================
# 41. LEARNED NETWORK VS FALLBACK
# =============================================================================

network_status <-
  feature_summary %>%
  
  dplyr::select(
    .data$season,
    .data$origin,
    .data$learned_network,
    .data$estimable_edges,
    .data$positive_edges
  )


all_forecasts <-
  all_forecasts %>%
  
  dplyr::left_join(
    
    network_status,
    
    by =
      c(
        "season",
        "origin"
      )
  )



network_performance <-
  all_forecasts %>%
  
  dplyr::group_by(
    .data$learned_network
  ) %>%
  
  dplyr::summarise(
    
    origins =
      dplyr::n_distinct(
        paste(
          .data$season,
          .data$origin
        )
      ),
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error,
          na.rm = TRUE
        )
      ),
    
    mean_WIS =
      mean(
        .data$WIS,
        na.rm = TRUE
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS,
        na.rm = TRUE
      ),
    
    coverage90 =
      mean(
        .data$covered90,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )



# =============================================================================
# 42. INTERACTION SUMMARY
# =============================================================================

interaction_summary <-
  parameter_summary %>%
  
  dplyr::filter(
    .data$variable ==
      "beta_interaction"
  ) %>%
  
  dplyr::mutate(
    
    positive_median =
      .data$median >
      0,
    
    cri90_above_zero =
      .data$q05 >
      0
  )



interaction_overall <-
  interaction_summary %>%
  
  dplyr::summarise(
    
    n_origins =
      dplyr::n(),
    
    mean_posterior_median =
      mean(
        .data$median
      ),
    
    median_posterior_median =
      stats::median(
        .data$median
      ),
    
    minimum_posterior_median =
      min(
        .data$median
      ),
    
    maximum_posterior_median =
      max(
        .data$median
      ),
    
    proportion_positive_median =
      mean(
        .data$positive_median
      ),
    
    proportion_90_CrI_above_zero =
      mean(
        .data$cri90_above_zero
      ),
    
    mean_probability_positive =
      mean(
        .data$probability_positive
      )
  )



# =============================================================================
# 43. ALL CORE PARAMETERS ACROSS ORIGINS
# =============================================================================

core_parameter_summary <-
  parameter_summary %>%
  
  dplyr::filter(
    
    .data$variable %in%
      
      c(
        "beta_memory",
        "beta_national",
        "beta_borrow",
        "beta_phase",
        "beta_interaction"
      )
  ) %>%
  
  dplyr::group_by(
    .data$variable
  ) %>%
  
  dplyr::summarise(
    
    n_origins =
      dplyr::n(),
    
    mean_median =
      mean(
        .data$median
      ),
    
    median_median =
      stats::median(
        .data$median
      ),
    
    min_median =
      min(
        .data$median
      ),
    
    max_median =
      max(
        .data$median
      ),
    
    proportion_positive =
      mean(
        .data$median >
          0
      ),
    
    proportion_90_CrI_above_zero =
      mean(
        .data$q05 >
          0
      ),
    
    proportion_90_CrI_below_zero =
      mean(
        .data$q95 <
          0
      ),
    
    .groups =
      "drop"
  )



# =============================================================================
# 44. CONVERGENCE OVERALL
# =============================================================================

convergence_overall <-
  convergence_summary %>%
  
  dplyr::summarise(
    
    total_fits =
      dplyr::n(),
    
    fits_passing =
      sum(
        .data$convergence_ok
      ),
    
    proportion_passing =
      mean(
        .data$convergence_ok
      ),
    
    maximum_Rhat =
      max(
        .data$max_rhat,
        na.rm = TRUE
      ),
    
    minimum_bulk_ESS =
      min(
        .data$min_ess_bulk,
        na.rm = TRUE
      ),
    
    minimum_tail_ESS =
      min(
        .data$min_ess_tail,
        na.rm = TRUE
      ),
    
    total_divergences =
      sum(
        .data$divergences,
        na.rm = TRUE
      ),
    
    total_treedepth_hits =
      sum(
        .data$treedepth_hits,
        na.rm = TRUE
      )
  )



# =============================================================================
# 45. SAVE FINAL TABLES
# =============================================================================

readr::write_csv(
  
  feature_summary,
  
  file.path(
    OUTPUT_DIR,
    "01_FEATURE_SUMMARY.csv"
  )
)


readr::write_csv(
  
  pairing_check,
  
  file.path(
    OUTPUT_DIR,
    "02_COMPLETENESS_CHECK.csv"
  )
)


readr::write_csv(
  
  convergence_summary,
  
  file.path(
    OUTPUT_DIR,
    "03_CONVERGENCE_BY_ORIGIN.csv"
  )
)


readr::write_csv(
  
  convergence_overall,
  
  file.path(
    OUTPUT_DIR,
    "04_CONVERGENCE_OVERALL.csv"
  )
)


readr::write_csv(
  
  parameter_summary,
  
  file.path(
    OUTPUT_DIR,
    "05_PARAMETER_SUMMARY_BY_ORIGIN.csv"
  )
)


readr::write_csv(
  
  interaction_summary,
  
  file.path(
    OUTPUT_DIR,
    "06_INTERACTION_BY_ORIGIN.csv"
  )
)


readr::write_csv(
  
  interaction_overall,
  
  file.path(
    OUTPUT_DIR,
    "07_INTERACTION_OVERALL.csv"
  )
)


readr::write_csv(
  
  core_parameter_summary,
  
  file.path(
    OUTPUT_DIR,
    "08_CORE_PARAMETER_SUMMARY.csv"
  )
)


readr::write_csv(
  
  origin_performance,
  
  file.path(
    OUTPUT_DIR,
    "09_PERFORMANCE_BY_ORIGIN.csv"
  )
)


readr::write_csv(
  
  overall_performance,
  
  file.path(
    OUTPUT_DIR,
    "10_OVERALL_PERFORMANCE.csv"
  )
)


readr::write_csv(
  
  season_performance,
  
  file.path(
    OUTPUT_DIR,
    "11_PERFORMANCE_BY_SEASON.csv"
  )
)


readr::write_csv(
  
  calibration_summary,
  
  file.path(
    OUTPUT_DIR,
    "12_CALIBRATION_OVERALL.csv"
  )
)


readr::write_csv(
  
  interval_miss_direction,
  
  file.path(
    OUTPUT_DIR,
    "13_INTERVAL_MISS_DIRECTION.csv"
  )
)


readr::write_csv(
  
  intensity_performance,
  
  file.path(
    OUTPUT_DIR,
    "14_PERFORMANCE_BY_CAUSAL_INTENSITY.csv"
  )
)


readr::write_csv(
  
  network_performance,
  
  file.path(
    OUTPUT_DIR,
    "15_PERFORMANCE_BY_NETWORK_STATUS.csv"
  )
)


readr::write_csv(
  
  all_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "16_ALL_STATE_FORECASTS.csv"
  )
)



# =============================================================================
# 46. CALIBRATION PLOT
# =============================================================================

p_calibration <-
  ggplot2::ggplot(
    
    calibration_summary,
    
    ggplot2::aes(
      x =
        nominal,
      y =
        observed
    )
  ) +
  
  ggplot2::geom_abline(
    slope =
      1,
    intercept =
      0,
    linetype =
      2
  ) +
  
  ggplot2::geom_point(
    size =
      3
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        interval
    ),
    
    nudge_y =
      0.025
  ) +
  
  ggplot2::coord_equal(
    
    xlim =
      c(
        0.4,
        1
      ),
    
    ylim =
      c(
        0.4,
        1
      )
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU M4 Predictive Interval Calibration",
    
    x =
      "Nominal coverage",
    
    y =
      "Observed coverage"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )


ggplot2::ggsave(
  
  filename =
    file.path(
      PLOT_DIR,
      "01_CALIBRATION.png"
    ),
  
  plot =
    p_calibration,
  
  width =
    7,
  
  height =
    6,
  
  dpi =
    300
)



# =============================================================================
# 47. WIS BY ORIGIN
# =============================================================================

p_wis <-
  ggplot2::ggplot(
    
    origin_performance,
    
    ggplot2::aes(
      
      x =
        origin,
      
      y =
        mean_WIS,
      
      group =
        season,
      
      linetype =
        season
    )
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU M4 WIS Across Forecast Origins",
    
    x =
      "Forecast origin",
    
    y =
      "Mean WIS",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )


ggplot2::ggsave(
  
  filename =
    file.path(
      PLOT_DIR,
      "02_WIS_BY_ORIGIN.png"
    ),
  
  plot =
    p_wis,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    300
)



# =============================================================================
# 48. INTERACTION BY ORIGIN
# =============================================================================

p_interaction <-
  ggplot2::ggplot(
    
    interaction_summary,
    
    ggplot2::aes(
      
      x =
        origin,
      
      y =
        median,
      
      ymin =
        q05,
      
      ymax =
        q95,
      
      group =
        season,
      
      linetype =
        season
    )
  ) +
  
  ggplot2::geom_hline(
    
    yintercept =
      0,
    
    linetype =
      2
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::geom_errorbar(
    width =
      0.5
  ) +
  
  ggplot2::labs(
    
    title =
      "Phase × Connectivity Interaction Across Origins",
    
    x =
      "Forecast origin",
    
    y =
      "Posterior median and 90% credible interval",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )


ggplot2::ggsave(
  
  filename =
    file.path(
      PLOT_DIR,
      "03_INTERACTION_BY_ORIGIN.png"
    ),
  
  plot =
    p_interaction,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    300
)



# =============================================================================
# 49. PERFORMANCE BY CAUSAL INTENSITY
# =============================================================================

p_intensity <-
  ggplot2::ggplot(
    
    intensity_performance,
    
    ggplot2::aes(
      
      x =
        intensity_group,
      
      y =
        coverage90,
      
      group =
        1
    )
  ) +
  
  ggplot2::geom_hline(
    
    yintercept =
      0.90,
    
    linetype =
      2
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point(
    size =
      3
  ) +
  
  ggplot2::labs(
    
    title =
      "M4 90% Coverage by Causal Epidemic Intensity",
    
    x =
      "Origin intensity group",
    
    y =
      "Observed 90% coverage"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )


ggplot2::ggsave(
  
  filename =
    file.path(
      PLOT_DIR,
      "04_COVERAGE_BY_CAUSAL_INTENSITY.png"
    ),
  
  plot =
    p_intensity,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    300
)



# =============================================================================
# 50. SESSION INFORMATION
# =============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      OUTPUT_DIR,
      "17_SESSION_INFO.txt"
    )
)



# =============================================================================
# 51. FINAL CONSOLE REPORT
# =============================================================================

cat(
  "\n\n============================================================\n"
)


cat(
  "PACE-FLU FINAL M4 VALIDATION COMPLETE\n"
)


cat(
  "============================================================\n"
)


cat(
  "\nCOMPLETENESS\n"
)


safe_print(
  pairing_check
)


cat(
  "\nCONVERGENCE OVERALL\n"
)


safe_print(
  convergence_overall
)


cat(
  "\nOVERALL PERFORMANCE\n"
)


safe_print(
  overall_performance
)


cat(
  "\nPERFORMANCE BY SEASON\n"
)


safe_print(
  season_performance
)


cat(
  "\nCALIBRATION\n"
)


safe_print(
  calibration_summary
)


cat(
  "\nINTERVAL MISS DIRECTION\n"
)


safe_print(
  interval_miss_direction
)


cat(
  "\nINTERACTION OVERALL\n"
)


safe_print(
  interaction_overall
)


cat(
  "\nCORE PARAMETER SUMMARY\n"
)


safe_print(
  core_parameter_summary
)


cat(
  "\nPERFORMANCE BY CAUSAL INTENSITY\n"
)


safe_print(
  intensity_performance
)


cat(
  "\nNETWORK STATUS PERFORMANCE\n"
)


safe_print(
  network_performance
)


cat(
  "\nOUTPUT DIRECTORY:\n"
)


cat(
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


cat(
  "\n============================================================\n"
)


cat(
  "FINAL M4 VALIDATION FINISHED.\n"
)


cat(
  "DO NOT MODIFY THE MODEL AFTER THIS RUN WITHOUT DEFINING\n"
)


cat(
  "THE CHANGE AS A NEW MODEL AND REVALIDATING IT.\n"
)


cat(
  "============================================================\n"
)


# =============================================================================
# END OF FINAL PACE-FLU M4 VALIDATION
# =============================================================================