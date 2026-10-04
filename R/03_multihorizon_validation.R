# =============================================================================
# PACE-FLU: MULTI-ORIGIN, MULTI-HORIZON RECURSIVE VALIDATION
# =============================================================================
#
# PURPOSE
# -------
# Determine whether the long-horizon recursive drift observed in the
# 2025/26 week-40 experiment is:
#
#   1. systematic across seasons;
#   2. dependent on forecast origin;
#   3. dependent on forecast horizon;
#   4. associated with particular epidemic conditions.
#
# FROZEN MODEL
# ------------
# PACE-FLU M4 is NOT modified.
#
# Origins:
#   20, 24, 28, 32, 36, 40
#
# Horizons:
#   1--12 weeks, subject to season end.
#
# Seasons:
#   2023/2024
#   2024/2025
#   2025/2026
#
# CRITICAL LEAKAGE RULE
# ---------------------
# At origin t:
#
#   previous seasons: all outcome-bearing rows
#   target season:    predictor rows <= t-1
#
# Thus outcomes are observed only through t.
#
# Week t initializes recursive prediction of t+1.
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
# 1. SETTINGS
# =============================================================================

SEED <- 20260926

set.seed(SEED)


STATE_DATA_FILE <-
  "data/nigeria_flu_weekly_by_state.csv"


OUTPUT_DIR <-
  "PACE_FLU_MULTI_ORIGIN_RECURSIVE_VALIDATION"


FIT_DIR <-
  file.path(
    OUTPUT_DIR,
    "fits"
  )


FIGURE_DIR <-
  file.path(
    OUTPUT_DIR,
    "figures"
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
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  STAN_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


SEASON_LEVELS <-
  c(
    "2023/2024",
    "2024/2025",
    "2025/2026"
  )


ORIGINS <-
  c(
    20L,
    24L,
    28L,
    32L,
    36L,
    40L
  )


MAX_HORIZON <- 12L


# Memory
L_MEMORY <- 8L
MEMORY_DECAY <- 0.35


MEMORY_WEIGHTS <-
  exp(
    -MEMORY_DECAY *
      0:(L_MEMORY - 1L)
  )


MEMORY_WEIGHTS <-
  MEMORY_WEIGHTS /
  sum(MEMORY_WEIGHTS)


MIN_EDGE_N <- 20L


# Stan
N_CHAINS <- 4L
N_PARALLEL_CHAINS <- 4L

N_WARMUP <- 1000L
N_SAMPLING <- 1000L

ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15L

STAN_INIT <- 0
INITIAL_STEP_SIZE <- 0.01


# Recursive posterior trajectories.
#
# 1000 is reasonable for this large diagnostic experiment.
# Increase to 2000 for the final version if desired.
N_RECURSIVE_DRAWS <- 1000L


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

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


# =============================================================================
# 3. INTERVAL SCORE / WIS
# =============================================================================

interval_score <- function(
    lower,
    upper,
    y,
    alpha
) {
  
  (upper - lower) +
    
    (2 / alpha) *
    (lower - y) *
    as.numeric(y < lower) +
    
    (2 / alpha) *
    (y - upper) *
    as.numeric(y > upper)
}


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
  
  ae <-
    abs(
      y - q50
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
    0.5 * ae +
      0.25 * IS50 +
      0.10 * IS80 +
      0.05 * IS90
  ) / 3.5
}


# =============================================================================
# 4. EMPIRICAL CRPS
# =============================================================================

empirical_crps <- function(
    draws,
    y
) {
  
  draws <-
    draws[
      is.finite(draws)
    ]
  
  
  draws <-
    sort(draws)
  
  
  D <-
    length(draws)
  
  
  if (
    D < 2L
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  term1 <-
    mean(
      abs(
        draws - y
      )
    )
  
  
  idx <-
    seq_len(D)
  
  
  term2 <-
    sum(
      (
        2 * idx -
          D -
          1
      ) *
        draws
    ) /
    D^2
  
  
  term1 - term2
}


# =============================================================================
# 5. DRAW SUMMARY
# =============================================================================

summarise_draws <- function(
    draws,
    actual
) {
  
  q <-
    stats::quantile(
      
      draws,
      
      probs =
        c(
          0.05,
          0.10,
          0.25,
          0.50,
          0.75,
          0.90,
          0.95
        ),
      
      names = FALSE,
      
      na.rm = TRUE,
      
      type = 8
    )
  
  
  tibble::tibble(
    
    q05 = q[1],
    q10 = q[2],
    q25 = q[3],
    q50 = q[4],
    q75 = q[5],
    q90 = q[6],
    q95 = q[7],
    
    actual =
      actual,
    
    error =
      q[4] -
      actual,
    
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
        draws,
        actual
      ),
    
    covered50 =
      actual >= q[3] &
      actual <= q[5],
    
    covered80 =
      actual >= q[2] &
      actual <= q[6],
    
    covered90 =
      actual >= q[1] &
      actual <= q[7],
    
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


# =============================================================================
# 6. READ DATA
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
# 7. BASIC DATA PREPARATION
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
    
    !is.na(.data$season),
    
    !is.na(.data$state),
    
    !is.na(.data$week),
    
    !is.na(.data$population),
    
    .data$population > 0,
    
    !is.na(.data$cases)
  ) %>%
  
  dplyr::arrange(
    .data$season_index,
    .data$state,
    .data$week
  )


# =============================================================================
# 8. STATE INDEX
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


cat(
  "\nNumber of states/FCT: ",
  N_STATES,
  "\n",
  sep = ""
)


if (
  N_STATES != 37L
) {
  
  warning(
    paste0(
      "Expected 37 states/FCT; found ",
      N_STATES
    )
  )
}


dat <-
  dat %>%
  
  dplyr::mutate(
    
    state_id =
      match(
        .data$state,
        STATE_LEVELS
      )
  )


# =============================================================================
# 9. GEOGRAPHY
# =============================================================================

dat <-
  dat %>%
  
  dplyr::mutate(
    
    zone =
      toupper(
        trimws(
          as.character(
            .data$zone
          )
        )
      ),
    
    macro_region =
      dplyr::case_when(
        
        .data$zone %in%
          c(
            "NW",
            "NE",
            "NC"
          ) ~
          "North",
        
        .data$zone %in%
          c(
            "SW",
            "SE",
            "SS"
          ) ~
          "South",
        
        TRUE ~
          NA_character_
      )
  )


if (
  any(
    is.na(
      dat$macro_region
    )
  )
) {
  
  stop(
    "Some states could not be assigned to North/South."
  )
}


state_lookup <-
  dat %>%
  
  dplyr::distinct(
    .data$state_id,
    .data$state,
    .data$zone,
    .data$macro_region,
    .data$population
  ) %>%
  
  dplyr::arrange(
    .data$state_id
  )


population_vector <-
  state_lookup$population


log_offset_vector <-
  log(
    population_vector /
      100000
  )


# =============================================================================
# 10. INCIDENCE
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
# 11. NATIONAL SIGNAL
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
# 12. MEMORY FEATURES
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
# 13. GROWTH AND NEXT-WEEK OUTCOME
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
# 14. CONNECTIVITY FUNCTION
# =============================================================================

estimate_connectivity <- function(
    train_data,
    N_STATES
) {
  
  edge_results <-
    list()
  
  
  counter <- 1L
  
  
  for (
    target_idx in
    seq_len(
      N_STATES
    )
  ) {
    
    target_data <-
      train_data %>%
      
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
        train_data %>%
        
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
          
          coefs <-
            stats::coef(
              edge_fit
            )
          
          
          if (
            "source_current" %in%
            names(coefs)
          ) {
            
            source_coefficient <-
              unname(
                coefs[
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
        counter + 1L
    }
  }
  
  
  edges <-
    dplyr::bind_rows(
      edge_results
    )
  
  
  W <-
    matrix(
      
      0,
      
      nrow =
        N_STATES,
      
      ncol =
        N_STATES
    )
  
  
  fallback_targets <- 0L
  
  
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
      sum(weights) >
      0
    ) {
      
      weights <-
        weights /
        sum(weights)
      
    } else {
      
      fallback_targets <-
        fallback_targets + 1L
      
      
      weights[] <-
        1 /
        (N_STATES - 1L)
    }
    
    
    W[
      target_idx,
      target_edges$source_id
    ] <-
      weights
  }
  
  
  diag(W) <- 0
  
  
  list(
    
    W =
      W,
    
    edges =
      edges,
    
    estimable_edges =
      sum(
        is.finite(
          edges$coefficient
        )
      ),
    
    positive_edges =
      sum(
        is.finite(
          edges$coefficient
        ) &
          edges$coefficient >
          0
      ),
    
    fallback_targets =
      fallback_targets
  )
}


# =============================================================================
# 15. BORROWING FUNCTION
# =============================================================================

calculate_borrowing <- function(
    data_rows,
    W,
    N_STATES
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
      nrow(block) !=
      N_STATES
    ) {
      
      next
    }
    
    
    borrowed <-
      as.numeric(
        
        W %*%
          block$log_incidence
      )
    
    
    result[
      indices
    ] <-
      borrowed[
        data_rows$state_id[
          indices
        ]
      ]
  }
  
  
  result
}


# =============================================================================
# 16. FROZEN M4 STAN MODEL
# =============================================================================

stan_code <- '

data {

  int<lower=1> N;

  int<lower=1> J;

  array[N]
    int<lower=0>
    y;

  array[N]
    int<lower=1,upper=J>
    state_id;

  vector[N] log_offset;

  vector[N] memory;

  vector[N] national;

  vector[N] borrowing;

  vector[N] phase;

  vector[N] borrow_phase;
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

  vector[N]
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
'


STAN_FILE <-
  file.path(
    STAN_DIR,
    "PACE_FLU_M4_MULTI_ORIGIN.stan"
  )


writeLines(
  stan_code,
  STAN_FILE
)


cat(
  "\nCompiling frozen M4...\n"
)


model <-
  cmdstanr::cmdstan_model(
    STAN_FILE
  )


# =============================================================================
# 17. STORAGE
# =============================================================================

all_forecasts <- list()

all_diagnostics <- list()

all_origin_metadata <- list()

all_parameter_summaries <- list()

all_national_forecasts <- list()

all_macro_forecasts <- list()


forecast_counter <- 1L
diagnostic_counter <- 1L
metadata_counter <- 1L
parameter_counter <- 1L
national_counter <- 1L
macro_counter <- 1L


# =============================================================================
# 18. LOOP OVER SEASONS AND ORIGINS
# =============================================================================

for (
  target_season in
  SEASON_LEVELS
) {
  
  target_season_index <-
    match(
      target_season,
      SEASON_LEVELS
    )
  
  
  cat(
    "\n\n============================================================\n"
  )
  
  cat(
    "SEASON: ",
    target_season,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  
  for (
    origin in
    ORIGINS
  ) {
    
    cat(
      "\n------------------------------------------------------------\n"
    )
    
    cat(
      "ORIGIN WEEK: ",
      origin,
      "\n",
      sep = ""
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    
    # =========================================================================
    # 18.1 FORECAST HORIZON
    # =========================================================================
    
    max_available_week <-
      max(
        dat$week[
          dat$season ==
            target_season
        ],
        na.rm = TRUE
      )
    
    
    horizon <-
      min(
        MAX_HORIZON,
        max_available_week -
          origin
      )
    
    
    if (
      horizon < 1L
    ) {
      
      next
    }
    
    
    forecast_weeks <-
      (
        origin + 1L
      ):
      (
        origin + horizon
      )
    
    
    # =========================================================================
    # 18.2 CHECKPOINT
    # =========================================================================
    
    season_tag <-
      gsub(
        "/",
        "_",
        target_season
      )
    
    
    checkpoint_file <-
      file.path(
        
        FIT_DIR,
        
        paste0(
          "recursive_",
          season_tag,
          "_origin_",
          origin,
          ".rds"
        )
      )
    
    
    # =========================================================================
    # 18.3 LEAKAGE-FREE TRAINING DATA
    # =========================================================================
    
    train <-
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
              origin - 1L
          )
      )
    
    
    # =========================================================================
    # CRITICAL EARLY-SEASON CHECK
    #
    # For 2023/24 there are no prior seasons in the challenge dataset.
    # We therefore need enough within-season observations to estimate M4.
    # =========================================================================
    
    if (
      nrow(train) <
      N_STATES * 10L
    ) {
      
      warning(
        paste0(
          "Insufficient training data for ",
          target_season,
          " origin ",
          origin,
          ". Skipping."
        )
      )
      
      next
    }
    
    
    # =========================================================================
    # 18.4 TRAIN-ONLY GROWTH SCALE
    # =========================================================================
    
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
    
    
    # =========================================================================
    # 18.5 CONNECTIVITY
    # =========================================================================
    
    connectivity <-
      estimate_connectivity(
        
        train,
        
        N_STATES
      )
    
    
    W <-
      connectivity$W
    
    
    # =========================================================================
    # 18.6 TRAINING BORROWING
    # =========================================================================
    
    train$borrowing <-
      calculate_borrowing(
        
        train,
        
        W,
        
        N_STATES
      )
    
    
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
      ) %>%
      
      dplyr::filter(
        
        stats::complete.cases(
          
          .data$cases_next,
          
          .data$memory,
          
          .data$national_log_incidence,
          
          .data$borrowing,
          
          .data$p_growth,
          
          .data$borrow_phase,
          
          .data$log_offset
        )
      )
    
    
    # =========================================================================
    # 18.7 STAN DATA
    # =========================================================================
    
    stan_data <-
      list(
        
        N =
          nrow(train),
        
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
        
        log_offset =
          as.vector(
            train$log_offset
          ),
        
        memory =
          as.vector(
            train$memory
          ),
        
        national =
          as.vector(
            train$national_log_incidence
          ),
        
        borrowing =
          as.vector(
            train$borrowing
          ),
        
        phase =
          as.vector(
            train$p_growth
          ),
        
        borrow_phase =
          as.vector(
            train$borrow_phase
          )
      )
    
    
    # =========================================================================
    # 18.8 FIT OR LOAD CHECKPOINT
    # =========================================================================
    
    if (
      file.exists(
        checkpoint_file
      )
    ) {
      
      cat(
        "Loading checkpoint...\n"
      )
      
      
      checkpoint <-
        readRDS(
          checkpoint_file
        )
      
      
      posterior_matrix <-
        checkpoint$posterior_matrix
      
      diagnostic_row <-
        checkpoint$diagnostic
      
      parameter_summary <-
        checkpoint$parameter_summary
      
      
    } else {
      
      cat(
        "Fitting M4...\n"
      )
      
      
      fit <-
        model$sample(
          
          data =
            stan_data,
          
          seed =
            SEED +
            target_season_index *
            1000L +
            origin,
          
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
      
      
      # =======================================================================
      # Diagnostics
      # =======================================================================
      
      sm <-
        fit$summary(
          
          variables =
            c(
              "alpha",
              "state_raw",
              "sigma_state",
              "beta_memory",
              "beta_national",
              "beta_borrow",
              "beta_phase",
              "beta_interaction",
              "log_phi"
            )
        )
      
      
      diagnostic_raw <-
        fit$diagnostic_summary()
      
      
      diagnostic_row <-
        tibble::tibble(
          
          season =
            target_season,
          
          origin =
            origin,
          
          max_rhat =
            max(
              sm$rhat,
              na.rm = TRUE
            ),
          
          min_bulk_ESS =
            min(
              sm$ess_bulk,
              na.rm = TRUE
            ),
          
          min_tail_ESS =
            min(
              sm$ess_tail,
              na.rm = TRUE
            ),
          
          divergences =
            sum(
              diagnostic_raw$num_divergent,
              na.rm = TRUE
            ),
          
          treedepth_hits =
            sum(
              diagnostic_raw$num_max_treedepth,
              na.rm = TRUE
            )
        )
      
      
      # =======================================================================
      # Parameter summary
      # =======================================================================
      
      # =======================================================================
      # Parameter summary
      #
      # Do NOT use probs= inside fit$summary().
      # Some cmdstanr/posterior versions pass it incorrectly to
      # posterior::summarise_draws().
      # =======================================================================
      
      parameter_variables <-
        c(
          "beta_memory",
          "beta_national",
          "beta_borrow",
          "beta_phase",
          "beta_interaction",
          "sigma_state",
          "log_phi"
        )
      
      
      parameter_draws <-
        fit$draws(
          variables = parameter_variables,
          format = "matrix"
        )
      
      
      parameter_summary <-
        lapply(
          
          parameter_variables,
          
          function(v) {
            
            x <-
              parameter_draws[, v]
            
            
            tibble::tibble(
              
              variable = v,
              
              mean =
                mean(x),
              
              sd =
                stats::sd(x),
              
              q05 =
                unname(
                  stats::quantile(
                    x,
                    0.05,
                    type = 8
                  )
                ),
              
              q50 =
                unname(
                  stats::quantile(
                    x,
                    0.50,
                    type = 8
                  )
                ),
              
              q95 =
                unname(
                  stats::quantile(
                    x,
                    0.95,
                    type = 8
                  )
                ),
              
              prob_positive =
                mean(
                  x > 0
                )
            )
          }
        ) %>%
        
        dplyr::bind_rows() %>%
        
        dplyr::mutate(
          
          season =
            target_season,
          
          origin =
            origin
        )
      
      
      # =======================================================================
      # Posterior matrix
      # =======================================================================
      
      posterior_matrix <-
        fit$draws(
          
          variables =
            c(
              "alpha",
              "alpha_state",
              "beta_memory",
              "beta_national",
              "beta_borrow",
              "beta_phase",
              "beta_interaction",
              "log_phi"
            ),
          
          format =
            "matrix"
        )
      
      
      saveRDS(
        
        list(
          
          posterior_matrix =
            posterior_matrix,
          
          diagnostic =
            diagnostic_row,
          
          parameter_summary =
            parameter_summary
        ),
        
        checkpoint_file
      )
      
      
      rm(fit)
      
      gc()
    }
    
    
    # =========================================================================
    # 18.9 SAVE DIAGNOSTICS
    # =========================================================================
    
    all_diagnostics[[
      diagnostic_counter
    ]] <-
      diagnostic_row
    
    
    diagnostic_counter <-
      diagnostic_counter + 1L
    
    
    all_parameter_summaries[[
      parameter_counter
    ]] <-
      parameter_summary
    
    
    parameter_counter <-
      parameter_counter + 1L
    
    
    # =========================================================================
    # 18.10 SELECT POSTERIOR DRAWS
    # =========================================================================
    
    TOTAL_DRAWS <-
      nrow(
        posterior_matrix
      )
    
    
    n_use <-
      min(
        N_RECURSIVE_DRAWS,
        TOTAL_DRAWS
      )
    
    
    draw_indices <-
      unique(
        round(
          seq(
            1,
            TOTAL_DRAWS,
            length.out =
              n_use
          )
        )
      )
    
    
    posterior_use <-
      posterior_matrix[
        draw_indices,
        ,
        drop = FALSE
      ]
    
    
    N_DRAWS <-
      nrow(
        posterior_use
      )
    
    
    # =========================================================================
    # 18.11 PARAMETERS
    # =========================================================================
    
    alpha_draw <-
      posterior_use[
        ,
        "alpha"
      ]
    
    
    beta_memory_draw <-
      posterior_use[
        ,
        "beta_memory"
      ]
    
    
    beta_national_draw <-
      posterior_use[
        ,
        "beta_national"
      ]
    
    
    beta_borrow_draw <-
      posterior_use[
        ,
        "beta_borrow"
      ]
    
    
    beta_phase_draw <-
      posterior_use[
        ,
        "beta_phase"
      ]
    
    
    beta_interaction_draw <-
      posterior_use[
        ,
        "beta_interaction"
      ]
    
    
    phi_draw <-
      0.05 +
      exp(
        posterior_use[
          ,
          "log_phi"
        ]
      )
    
    
    alpha_state_draw <-
      matrix(
        
        NA_real_,
        
        nrow =
          N_DRAWS,
        
        ncol =
          N_STATES
      )
    
    
    for (
      state_idx in
      seq_len(
        N_STATES
      )
    ) {
      
      alpha_state_draw[
        ,
        state_idx
      ] <-
        
        posterior_use[
          ,
          paste0(
            "alpha_state[",
            state_idx,
            "]"
          )
        ]
    }
    
    
    # =========================================================================
    # 18.12 OBSERVED TARGET-SEASON DATA
    # =========================================================================
    
    target_data <-
      dat %>%
      
      dplyr::filter(
        .data$season ==
          target_season
      ) %>%
      
      dplyr::arrange(
        .data$state_id,
        .data$week
      )
    
    
    # =========================================================================
    # 18.13 TRAJECTORY ARRAY
    # =========================================================================
    
    trajectory_length <-
      origin +
      horizon
    
    
    case_array <-
      array(
        
        NA_real_,
        
        dim =
          c(
            N_DRAWS,
            N_STATES,
            trajectory_length
          )
      )
    
    
    # =========================================================================
    # Fill only observations available at origin.
    # =========================================================================
    
    for (
      state_idx in
      seq_len(
        N_STATES
      )
    ) {
      
      observed_history <-
        target_data %>%
        
        dplyr::filter(
          
          .data$state_id ==
            state_idx,
          
          .data$week <=
            origin
        ) %>%
        
        dplyr::arrange(
          .data$week
        )
      
      
      if (
        nrow(
          observed_history
        ) !=
        origin
      ) {
        
        stop(
          paste0(
            "Incomplete history: ",
            target_season,
            ", origin ",
            origin,
            ", state ",
            state_idx
          )
        )
      }
      
      
      case_array[
        ,
        state_idx,
        1:origin
      ] <-
        
        matrix(
          
          observed_history$cases,
          
          nrow =
            N_DRAWS,
          
          ncol =
            origin,
          
          byrow =
            TRUE
        )
    }
    
    
    # =========================================================================
    # 18.14 RECURSIVE FORECASTING
    # =========================================================================
    
    for (
      future_week in
      forecast_weeks
    ) {
      
      cat(
        "  Forecasting week ",
        future_week,
        " (h=",
        future_week - origin,
        ")\n",
        sep = ""
      )
      
      
      previous_week <-
        future_week - 1L
      
      
      for (
        draw_idx in
        seq_len(
          N_DRAWS
        )
      ) {
        
        # ---------------------------------------------------------------------
        # Current and lagged log incidence
        # ---------------------------------------------------------------------
        
        current_cases <-
          case_array[
            draw_idx,
            ,
            previous_week
          ]
        
        
        previous_cases <-
          case_array[
            draw_idx,
            ,
            previous_week - 1L
          ]
        
        
        current_log_incidence <-
          log1p(
            
            100000 *
              current_cases /
              population_vector
          )
        
        
        previous_log_incidence <-
          log1p(
            
            100000 *
              previous_cases /
              population_vector
          )
        
        
        # ---------------------------------------------------------------------
        # Memory
        # ---------------------------------------------------------------------
        
        memory_current <-
          numeric(
            N_STATES
          )
        
        
        for (
          state_idx in
          seq_len(
            N_STATES
          )
        ) {
          
          recent_weeks <-
            previous_week -
            0:(L_MEMORY - 1L)
          
          
          recent_cases <-
            case_array[
              draw_idx,
              state_idx,
              recent_weeks
            ]
          
          
          recent_log_incidence <-
            log1p(
              
              100000 *
                recent_cases /
                population_vector[
                  state_idx
                ]
            )
          
          
          memory_current[
            state_idx
          ] <-
            sum(
              MEMORY_WEIGHTS *
                recent_log_incidence
            )
        }
        
        
        # ---------------------------------------------------------------------
        # Phase
        # ---------------------------------------------------------------------
        
        growth_current <-
          current_log_incidence -
          previous_log_incidence
        
        
        phase_current <-
          plogis(
            
            growth_current /
              growth_scale
          )
        
        
        # ---------------------------------------------------------------------
        # Connectivity borrowing
        # ---------------------------------------------------------------------
        
        borrowing_current <-
          as.numeric(
            
            W %*%
              current_log_incidence
          )
        
        
        borrow_phase_current <-
          borrowing_current *
          phase_current
        
        
        # ---------------------------------------------------------------------
        # National signal
        # ---------------------------------------------------------------------
        
        national_current <-
          log1p(
            
            100000 *
              sum(
                current_cases
              ) /
              sum(
                population_vector
              )
          )
        
        
        # ---------------------------------------------------------------------
        # M4 predictor
        # ---------------------------------------------------------------------
        
        eta <-
          
          log_offset_vector +
          
          alpha_draw[
            draw_idx
          ] +
          
          alpha_state_draw[
            draw_idx,
          ] +
          
          beta_memory_draw[
            draw_idx
          ] *
          memory_current +
          
          beta_national_draw[
            draw_idx
          ] *
          national_current +
          
          beta_borrow_draw[
            draw_idx
          ] *
          borrowing_current +
          
          beta_phase_draw[
            draw_idx
          ] *
          phase_current +
          
          beta_interaction_draw[
            draw_idx
          ] *
          borrow_phase_current
        
        
        if (
          any(
            !is.finite(
              eta
            )
          )
        ) {
          
          stop(
            paste0(
              "Non-finite eta: ",
              target_season,
              " origin ",
              origin,
              " week ",
              future_week
            )
          )
        }
        
        
        mu <-
          exp(
            eta
          )
        
        
        if (
          any(
            !is.finite(
              mu
            )
          )
        ) {
          
          stop(
            paste0(
              "Non-finite mu: ",
              target_season,
              " origin ",
              origin,
              " week ",
              future_week
            )
          )
        }
        
        
        case_array[
          draw_idx,
          ,
          future_week
        ] <-
          
          stats::rnbinom(
            
            n =
              N_STATES,
            
            size =
              phi_draw[
                draw_idx
              ],
            
            mu =
              mu
          )
      }
    }
    
    
    # =========================================================================
    # 18.15 SCORE STATE FORECASTS
    # =========================================================================
    
    for (
      future_week in
      forecast_weeks
    ) {
      
      h <-
        future_week -
        origin
      
      
      for (
        state_idx in
        seq_len(
          N_STATES
        )
      ) {
        
        actual_row <-
          target_data %>%
          
          dplyr::filter(
            
            .data$state_id ==
              state_idx,
            
            .data$week ==
              future_week
          )
        
        
        if (
          nrow(
            actual_row
          ) !=
          1L
        ) {
          
          next
        }
        
        
        forecast_draws <-
          case_array[
            ,
            state_idx,
            future_week
          ]
        
        
        forecast_summary <-
          summarise_draws(
            
            forecast_draws,
            
            actual_row$cases[
              1
            ]
          )
        
        
        all_forecasts[[
          forecast_counter
        ]] <-
          
          tibble::tibble(
            
            season =
              target_season,
            
            origin =
              origin,
            
            target_week =
              future_week,
            
            horizon =
              h,
            
            state_id =
              state_idx,
            
            state =
              state_lookup$state[
                state_idx
              ],
            
            zone =
              state_lookup$zone[
                state_idx
              ],
            
            macro_region =
              state_lookup$macro_region[
                state_idx
              ]
          ) %>%
          
          dplyr::bind_cols(
            forecast_summary
          )
        
        
        forecast_counter <-
          forecast_counter + 1L
      }
      
      
      # =======================================================================
      # 18.16 NATIONAL FORECAST
      # =======================================================================
      
      national_draws <-
        rowSums(
          
          case_array[
            ,
            ,
            future_week,
            drop = FALSE
          ][
            ,
            ,
            1
          ]
        )
      
      
      national_actual <-
        target_data %>%
        
        dplyr::filter(
          .data$week ==
            future_week
        ) %>%
        
        dplyr::summarise(
          
          cases =
            sum(
              .data$cases
            )
        ) %>%
        
        dplyr::pull(
          .data$cases
        )
      
      
      all_national_forecasts[[
        national_counter
      ]] <-
        
        tibble::tibble(
          
          season =
            target_season,
          
          origin =
            origin,
          
          target_week =
            future_week,
          
          horizon =
            h
        ) %>%
        
        dplyr::bind_cols(
          
          summarise_draws(
            
            national_draws,
            
            national_actual
          )
        )
      
      
      national_counter <-
        national_counter + 1L
      
      
      # =======================================================================
      # 18.17 NORTH/SOUTH
      # =======================================================================
      
      for (
        macro_name in
        c(
          "North",
          "South"
        )
      ) {
        
        state_indices <-
          which(
            state_lookup$macro_region ==
              macro_name
          )
        
        
        macro_draws <-
          rowSums(
            
            case_array[
              ,
              state_indices,
              future_week,
              drop = FALSE
            ][
              ,
              ,
              1
            ]
          )
        
        
        macro_actual <-
          target_data %>%
          
          dplyr::filter(
            
            .data$week ==
              future_week,
            
            .data$state_id %in%
              state_indices
          ) %>%
          
          dplyr::summarise(
            
            cases =
              sum(
                .data$cases
              )
          ) %>%
          
          dplyr::pull(
            .data$cases
          )
        
        
        all_macro_forecasts[[
          macro_counter
        ]] <-
          
          tibble::tibble(
            
            season =
              target_season,
            
            origin =
              origin,
            
            target_week =
              future_week,
            
            horizon =
              h,
            
            region =
              macro_name
          ) %>%
          
          dplyr::bind_cols(
            
            summarise_draws(
              
              macro_draws,
              
              macro_actual
            )
          )
        
        
        macro_counter <-
          macro_counter + 1L
      }
    }
    
    
    # =========================================================================
    # 18.18 ORIGIN METADATA
    # =========================================================================
    
    origin_observed <-
      target_data %>%
      
      dplyr::filter(
        .data$week ==
          origin
      )
    
    
    origin_national_cases <-
      sum(
        origin_observed$cases
      )
    
    
    origin_national_incidence <-
      100000 *
      origin_national_cases /
      sum(
        population_vector
      )
    
    
    all_origin_metadata[[
      metadata_counter
    ]] <-
      
      tibble::tibble(
        
        season =
          target_season,
        
        origin =
          origin,
        
        horizon =
          horizon,
        
        training_rows =
          nrow(
            train
          ),
        
        growth_scale =
          growth_scale,
        
        estimable_edges =
          connectivity$estimable_edges,
        
        positive_edges =
          connectivity$positive_edges,
        
        fallback_targets =
          connectivity$fallback_targets,
        
        origin_national_cases =
          origin_national_cases,
        
        origin_national_incidence =
          origin_national_incidence
      )
    
    
    metadata_counter <-
      metadata_counter + 1L
    
    
    # =========================================================================
    # MEMORY CLEANUP
    # =========================================================================
    
    rm(
      case_array,
      posterior_matrix,
      posterior_use
    )
    
    gc()
  }
}


# =============================================================================
# 19. COMBINE RESULTS
# =============================================================================

state_forecasts <-
  dplyr::bind_rows(
    all_forecasts
  )


diagnostics <-
  dplyr::bind_rows(
    all_diagnostics
  )


origin_metadata <-
  dplyr::bind_rows(
    all_origin_metadata
  )


parameter_summaries <-
  dplyr::bind_rows(
    all_parameter_summaries
  )


national_forecasts <-
  dplyr::bind_rows(
    all_national_forecasts
  )


macro_forecasts <-
  dplyr::bind_rows(
    all_macro_forecasts
  )


# =============================================================================
# 20. SAVE RAW RESULTS
# =============================================================================

readr::write_csv(
  
  state_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "01_ALL_STATE_RECURSIVE_FORECASTS.csv"
  )
)


readr::write_csv(
  
  national_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "02_ALL_NATIONAL_RECURSIVE_FORECASTS.csv"
  )
)


readr::write_csv(
  
  macro_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "03_ALL_NORTH_SOUTH_RECURSIVE_FORECASTS.csv"
  )
)


readr::write_csv(
  
  diagnostics,
  
  file.path(
    OUTPUT_DIR,
    "04_MODEL_DIAGNOSTICS.csv"
  )
)


readr::write_csv(
  
  origin_metadata,
  
  file.path(
    OUTPUT_DIR,
    "05_ORIGIN_METADATA.csv"
  )
)


readr::write_csv(
  
  parameter_summaries,
  
  file.path(
    OUTPUT_DIR,
    "06_PARAMETER_SUMMARIES.csv"
  )
)


# =============================================================================
# 21. PERFORMANCE BY HORIZON
# =============================================================================

performance_by_horizon <-
  state_forecasts %>%
  
  dplyr::group_by(
    .data$horizon
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    median_AE =
      stats::median(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    median_bias =
      stats::median(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    median_WIS =
      stats::median(
        .data$WIS
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      ),
    
    mean_width90 =
      mean(
        .data$width90
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  performance_by_horizon,
  
  file.path(
    OUTPUT_DIR,
    "07_PERFORMANCE_BY_HORIZON.csv"
  )
)


# =============================================================================
# 22. PERFORMANCE BY SEASON AND HORIZON
# =============================================================================

performance_season_horizon <-
  state_forecasts %>%
  
  dplyr::group_by(
    .data$season,
    .data$horizon
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  performance_season_horizon,
  
  file.path(
    OUTPUT_DIR,
    "08_PERFORMANCE_BY_SEASON_AND_HORIZON.csv"
  )
)


# =============================================================================
# 23. PERFORMANCE BY ORIGIN
# =============================================================================

performance_by_origin <-
  state_forecasts %>%
  
  dplyr::group_by(
    .data$season,
    .data$origin
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  performance_by_origin,
  
  file.path(
    OUTPUT_DIR,
    "09_PERFORMANCE_BY_ORIGIN.csv"
  )
)


# =============================================================================
# 24. NATIONAL PERFORMANCE BY HORIZON
# =============================================================================

national_performance <-
  national_forecasts %>%
  
  dplyr::group_by(
    .data$horizon
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  national_performance,
  
  file.path(
    OUTPUT_DIR,
    "10_NATIONAL_PERFORMANCE_BY_HORIZON.csv"
  )
)


# =============================================================================
# 25. NORTH/SOUTH PERFORMANCE
# =============================================================================

macro_performance <-
  macro_forecasts %>%
  
  dplyr::group_by(
    .data$region,
    .data$horizon
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  macro_performance,
  
  file.path(
    OUTPUT_DIR,
    "11_NORTH_SOUTH_PERFORMANCE_BY_HORIZON.csv"
  )
)


# =============================================================================
# 26. BIAS-DIRECTION DIAGNOSTIC
# =============================================================================

bias_direction <-
  state_forecasts %>%
  
  dplyr::group_by(
    .data$horizon
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    proportion_overprediction =
      mean(
        .data$error > 0
      ),
    
    proportion_underprediction =
      mean(
        .data$error < 0
      ),
    
    mean_bias =
      mean(
        .data$error
      ),
    
    median_bias =
      stats::median(
        .data$error
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  bias_direction,
  
  file.path(
    OUTPUT_DIR,
    "12_BIAS_DIRECTION_BY_HORIZON.csv"
  )
)


# =============================================================================
# 27. OVERALL SUMMARY
# =============================================================================

overall_summary <-
  state_forecasts %>%
  
  dplyr::summarise(
    
    n_forecasts =
      dplyr::n(),
    
    n_seasons =
      dplyr::n_distinct(
        .data$season
      ),
    
    n_origins =
      dplyr::n_distinct(
        paste(
          .data$season,
          .data$origin
        )
      ),
    
    max_horizon =
      max(
        .data$horizon
      ),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    median_AE =
      stats::median(
        .data$absolute_error
      ),
    
    RMSE =
      sqrt(
        mean(
          .data$squared_error
        )
      ),
    
    bias =
      mean(
        .data$error
      ),
    
    median_bias =
      stats::median(
        .data$error
      ),
    
    mean_WIS =
      mean(
        .data$WIS
      ),
    
    median_WIS =
      stats::median(
        .data$WIS
      ),
    
    mean_CRPS =
      mean(
        .data$CRPS
      ),
    
    coverage50 =
      mean(
        .data$covered50
      ),
    
    coverage80 =
      mean(
        .data$covered80
      ),
    
    coverage90 =
      mean(
        .data$covered90
      )
  )


readr::write_csv(
  
  overall_summary,
  
  file.path(
    OUTPUT_DIR,
    "13_OVERALL_RECURSIVE_SUMMARY.csv"
  )
)


# =============================================================================
# 28. DIAGNOSTIC FLAGS
#
# These are descriptive, not model-selection rules.
# =============================================================================

horizon_diagnostic <-
  performance_by_horizon %>%
  
  dplyr::mutate(
    
    relative_MAE =
      .data$MAE /
      dplyr::first(
        .data$MAE
      ),
    
    relative_WIS =
      .data$mean_WIS /
      dplyr::first(
        .data$mean_WIS
      ),
    
    absolute_calibration_error90 =
      abs(
        .data$coverage90 -
          0.90
      )
  )


readr::write_csv(
  
  horizon_diagnostic,
  
  file.path(
    OUTPUT_DIR,
    "14_HORIZON_DIAGNOSTIC.csv"
  )
)


# =============================================================================
# 29. FIGURE: MAE BY HORIZON
# =============================================================================

p_mae <-
  ggplot2::ggplot(
    
    performance_by_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = MAE
    )
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::scale_x_continuous(
    breaks =
      1:MAX_HORIZON
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU Recursive Forecast Error by Horizon",
    
    subtitle =
      "Pooled across seasons and forecast origins",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Mean absolute error"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "01_MAE_BY_HORIZON.png"
  ),
  
  p_mae,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 30. FIGURE: WIS BY HORIZON
# =============================================================================

p_wis <-
  ggplot2::ggplot(
    
    performance_by_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = mean_WIS
    )
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::scale_x_continuous(
    breaks =
      1:MAX_HORIZON
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU Probabilistic Accuracy by Horizon",
    
    subtitle =
      "Multi-origin recursive validation",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Mean WIS"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "02_WIS_BY_HORIZON.png"
  ),
  
  p_wis,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 31. FIGURE: BIAS BY HORIZON
# =============================================================================

p_bias <-
  ggplot2::ggplot(
    
    performance_by_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = bias
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::scale_x_continuous(
    breaks =
      1:MAX_HORIZON
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU Recursive Forecast Bias",
    
    subtitle =
      "Positive values indicate overprediction",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Mean forecast error"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "03_BIAS_BY_HORIZON.png"
  ),
  
  p_bias,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 32. FIGURE: SEASON-SPECIFIC BIAS
# =============================================================================

p_season_bias <-
  ggplot2::ggplot(
    
    performance_season_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = bias,
      group = season
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  
  ggplot2::geom_line(
    
    ggplot2::aes(
      linetype = season
    )
  ) +
  
  ggplot2::geom_point() +
  
  ggplot2::scale_x_continuous(
    breaks =
      1:MAX_HORIZON
  ) +
  
  ggplot2::labs(
    
    title =
      "Recursive Forecast Bias by Season",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Mean forecast error",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "04_BIAS_BY_SEASON_AND_HORIZON.png"
  ),
  
  p_season_bias,
  
  width = 9,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 33. FIGURE: 90% COVERAGE
# =============================================================================

p_coverage <-
  ggplot2::ggplot(
    
    performance_by_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = coverage90
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept = 0.90,
    linetype = 2
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::scale_x_continuous(
    breaks =
      1:MAX_HORIZON
  ) +
  
  ggplot2::scale_y_continuous(
    limits =
      c(
        0,
        1
      )
  ) +
  
  ggplot2::labs(
    
    title =
      "Empirical 90% Forecast-Interval Coverage",
    
    subtitle =
      "Dashed line represents nominal 90% coverage",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Empirical coverage"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "05_COVERAGE90_BY_HORIZON.png"
  ),
  
  p_coverage,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 34. SESSION INFORMATION
# =============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      OUTPUT_DIR,
      "15_SESSION_INFO.txt"
    )
)


# =============================================================================
# 35. PRINT RESULTS
# =============================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "MULTI-ORIGIN RECURSIVE VALIDATION COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nOVERALL SUMMARY\n"
)

print(
  overall_summary
)


cat(
  "\nPERFORMANCE BY HORIZON\n"
)

print(
  performance_by_horizon,
  n = Inf
)


cat(
  "\nBIAS DIRECTION\n"
)

print(
  bias_direction,
  n = Inf
)


cat(
  "\nMODEL DIAGNOSTICS\n"
)

print(
  diagnostics,
  n = Inf
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
  "END\n"
)

cat(
  "============================================================\n"
)