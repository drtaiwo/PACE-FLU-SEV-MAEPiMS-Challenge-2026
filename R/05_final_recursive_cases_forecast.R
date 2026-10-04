# =============================================================================
# PACE-FLU FINAL RECURSIVE CASE FORECASTING PIPELINE
#
# Frozen model: M4
# Season:       2025/2026
# Origin:       epidemiological week 40
# Forecast:     weeks 41--52
#
# IMPORTANT:
#   * Weeks 41--52 are NEVER used for fitting.
#   * Training rows in 2025/26 stop at predictor week 39.
#   * Week 40 is the observed forecast-origin state.
#   * Connectivity is learned at origin and frozen.
#   * Memory, phase, borrowing and national activity evolve recursively.
#   * State trajectories are aggregated DRAW-BY-DRAW.
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

SEED <- 20260926
set.seed(SEED)

STATE_DATA_FILE <-
  "data/nigeria_flu_weekly_by_state.csv"

METADATA_FILE <-
  "data/nigeria_flu_state_metadata.csv"

OUTPUT_DIR <-
  "PACE_FLU_FINAL_RECURSIVE_CASES"

STAN_DIR <-
  file.path(OUTPUT_DIR, "stan")

FIGURE_DIR <-
  file.path(OUTPUT_DIR, "figures")

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  STAN_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


TARGET_SEASON <- "2025/2026"
TARGET_SEASON_INDEX <- 3L

FORECAST_ORIGIN <- 40L

FORECAST_WEEKS <- 41:52

HORIZON <- length(FORECAST_WEEKS)


SEASON_LEVELS <- c(
  "2023/2024",
  "2024/2025",
  "2025/2026"
)


L_MEMORY <- 8L
MEMORY_DECAY <- 0.35

MEMORY_WEIGHTS <-
  exp(
    -MEMORY_DECAY *
      (0:(L_MEMORY - 1L))
  )

MEMORY_WEIGHTS <-
  MEMORY_WEIGHTS /
  sum(MEMORY_WEIGHTS)


MIN_EDGE_N <- 20L


# MCMC
N_CHAINS <- 4L
N_PARALLEL_CHAINS <- 4L

N_WARMUP <- 1000L
N_SAMPLING <- 1000L

ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15L

STAN_INIT <- 0
INITIAL_STEP_SIZE <- 0.01


# Number of posterior trajectories used recursively.
# 2000 is usually enough for stable quantiles while remaining manageable.
N_RECURSIVE_DRAWS <- 2000L


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

normalize_season <- function(x) {
  
  x <- as.character(x)
  
  dplyr::case_when(
    
    x %in% c(
      "2023/24",
      "2023/2024"
    ) ~ "2023/2024",
    
    x %in% c(
      "2024/25",
      "2024/2025"
    ) ~ "2024/2025",
    
    x %in% c(
      "2025/26",
      "2025/2026"
    ) ~ "2025/2026",
    
    TRUE ~ NA_character_
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


safe_print <- function(x) {
  
  print(
    as.data.frame(x),
    row.names = FALSE
  )
}


check_finite <- function(x, name) {
  
  if (
    any(!is.finite(x))
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
# 3. INTERVAL SCORE AND WIS
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
    abs(y - q50)
  
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

empirical_crps <- function(draws, y) {
  
  draws <-
    draws[
      is.finite(draws)
    ]
  
  draws <-
    sort(draws)
  
  D <-
    length(draws)
  
  if (D < 2L) {
    return(NA_real_)
  }
  
  term1 <-
    mean(
      abs(draws - y)
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
# 5. QUANTILE FUNCTION
# =============================================================================

summarise_draws <- function(
    draws,
    actual = NA_real_
) {
  
  q <-
    stats::quantile(
      
      draws,
      
      probs = c(
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95
      ),
      
      na.rm = TRUE,
      
      names = FALSE,
      
      type = 8
    )
  
  result <-
    tibble::tibble(
      
      q05 = q[1],
      q10 = q[2],
      q25 = q[3],
      q50 = q[4],
      q75 = q[5],
      q90 = q[6],
      q95 = q[7]
    )
  
  if (
    is.finite(actual)
  ) {
    
    result <-
      result %>%
      
      dplyr::mutate(
        
        actual = actual,
        
        error =
          .data$q50 -
          actual,
        
        absolute_error =
          abs(
            .data$q50 -
              actual
          ),
        
        squared_error =
          (
            .data$q50 -
              actual
          )^2,
        
        WIS =
          calculate_wis(
            
            actual,
            
            .data$q05,
            .data$q10,
            .data$q25,
            .data$q50,
            .data$q75,
            .data$q90,
            .data$q95
          ),
        
        CRPS =
          empirical_crps(
            draws,
            actual
          ),
        
        covered50 =
          actual >= .data$q25 &
          actual <= .data$q75,
        
        covered80 =
          actual >= .data$q10 &
          actual <= .data$q90,
        
        covered90 =
          actual >= .data$q05 &
          actual <= .data$q95,
        
        width50 =
          .data$q75 -
          .data$q25,
        
        width80 =
          .data$q90 -
          .data$q10,
        
        width90 =
          .data$q95 -
          .data$q05
      )
  }
  
  result
}


# =============================================================================
# 6. READ DATA
# =============================================================================

if (
  !file.exists(STATE_DATA_FILE)
) {
  stop(
    paste0(
      "Cannot find ",
      STATE_DATA_FILE
    )
  )
}


state_raw <-
  readr::read_csv(
    STATE_DATA_FILE,
    show_col_types = FALSE
  )


metadata <-
  readr::read_csv(
    METADATA_FILE,
    show_col_types = FALSE
  )


cat("\nSTATE DATA COLUMNS\n")
print(names(state_raw))

cat("\nMETADATA COLUMNS\n")
print(names(metadata))


# =============================================================================
# 7. PREPARE STATE DATA
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
    unique(dat$state)
  )


N_STATES <-
  length(STATE_LEVELS)


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
# 9. GEOGRAPHIC INFORMATION
#
# Prefer zone already present in state data.
# Metadata is used only if needed.
# =============================================================================

if (
  !"zone" %in% names(dat)
) {
  
  zone_candidates <-
    c(
      "zone",
      "geopolitical_zone",
      "geopolitical.zone",
      "region"
    )
  
  zone_name <-
    zone_candidates[
      zone_candidates %in%
        names(metadata)
    ][1]
  
  if (
    is.na(zone_name)
  ) {
    stop(
      "Could not identify geopolitical-zone column."
    )
  }
  
  metadata_small <-
    metadata %>%
    
    dplyr::transmute(
      
      state = .data$state,
      
      zone =
        .data[[zone_name]]
    )
  
  dat <-
    dat %>%
    
    dplyr::left_join(
      metadata_small,
      by = "state"
    )
}


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
          c("NW", "NE", "NC") ~
          "North",
        
        .data$zone %in%
          c("SW", "SE", "SS") ~
          "South",
        
        TRUE ~
          NA_character_
      )
  )


if (
  any(
    is.na(dat$macro_region)
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
# 11. NATIONAL EPIDEMIC SIGNAL
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
    
    .groups = "drop"
  )


dat <-
  dat %>%
  
  dplyr::left_join(
    
    national_data,
    
    by = c(
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
    seq_len(L_MEMORY)
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
# 13. NEXT-WEEK OUTCOME AND GROWTH
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
# 14. LEAKAGE-FREE TRAINING DATA
#
# Previous seasons:
# all outcome-bearing rows.
#
# 2025/26:
# predictor week <= 39
# so cases_next <= week 40.
# =============================================================================

train <-
  dat %>%
  
  dplyr::filter(
    
    !is.na(.data$cases_next),
    
    !is.na(.data$log_incidence_next),
    
    !is.na(.data$memory),
    
    !is.na(.data$growth_raw),
    
    (
      .data$season_index <
        TARGET_SEASON_INDEX
    ) |
      
      (
        .data$season_index ==
          TARGET_SEASON_INDEX &
          
          .data$week <=
          FORECAST_ORIGIN - 1L
      )
  )


# =============================================================================
# 15. OBSERVED WEEK-40 ORIGIN STATE
# =============================================================================

origin_data <-
  dat %>%
  
  dplyr::filter(
    
    .data$season ==
      TARGET_SEASON,
    
    .data$week ==
      FORECAST_ORIGIN
  ) %>%
  
  dplyr::arrange(
    .data$state_id
  )


if (
  nrow(origin_data) !=
  N_STATES
) {
  
  stop(
    "Week-40 origin does not contain all states."
  )
}


# =============================================================================
# 16. GROWTH SCALE
# =============================================================================

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


origin_data <-
  origin_data %>%
  
  dplyr::mutate(
    
    p_growth =
      plogis(
        .data$growth_raw /
          growth_scale
      )
  )


cat(
  "\nGrowth scale = ",
  growth_scale,
  "\n",
  sep = ""
)


# =============================================================================
# 17. CONNECTIVITY ESTIMATION
# =============================================================================

estimate_connectivity <- function(
    train_data
) {
  
  edge_results <- list()
  
  counter <- 1L
  
  
  for (
    target_idx in
    seq_len(N_STATES)
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
      seq_len(N_STATES)
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
          
          by = c(
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
        nrow(edge_data)
      
      
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
            
            silent = TRUE
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
      
      
      edge_results[[counter]] <-
        
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
      nrow = N_STATES,
      ncol = N_STATES
    )
  
  
  for (
    target_idx in
    seq_len(N_STATES)
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
  
  
  if (
    max(
      abs(
        rowSums(W) -
        1
      )
    ) >
    1e-8
  ) {
    
    stop(
      "Connectivity matrix does not sum to one."
    )
  }
  
  
  list(
    
    W = W,
    
    edges = edges,
    
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
          edges$coefficient > 0
      )
  )
}


connectivity <-
  estimate_connectivity(
    train
  )


W <-
  connectivity$W


cat(
  "\nEstimable edges: ",
  connectivity$estimable_edges,
  "\nPositive edges: ",
  connectivity$positive_edges,
  "\n",
  sep = ""
)


# =============================================================================
# 18. BORROWING FOR TRAINING DATA
# =============================================================================

calculate_borrowing <- function(
    data_rows,
    W
) {
  
  result <-
    rep(
      NA_real_,
      nrow(data_rows)
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
      nrow(blocks)
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
    
    
    source_signal <-
      block$log_incidence
    
    
    borrowed <-
      as.numeric(
        W %*%
          source_signal
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


train$borrowing <-
  calculate_borrowing(
    train,
    W
  )


origin_data$borrowing <-
  as.numeric(
    W %*%
      origin_data$log_incidence
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
  )


origin_data <-
  origin_data %>%
  
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


train <-
  train %>%
  
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


# =============================================================================
# 19. FINITE CHECKS
# =============================================================================

predictors <- c(
  "memory",
  "national_log_incidence",
  "borrowing",
  "p_growth",
  "borrow_phase",
  "log_offset"
)


for (
  variable_name in
  predictors
) {
  
  check_finite(
    train[[variable_name]],
    paste0(
      "training ",
      variable_name
    )
  )
  
  check_finite(
    origin_data[[variable_name]],
    paste0(
      "origin ",
      variable_name
    )
  )
}


# =============================================================================
# 20. STAN MODEL
#
# This is the frozen M4.
#
# No future y_rep is needed here.
# Recursive forecasting will be performed in R using posterior parameter draws.
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
    "PACE_FLU_FINAL_M4_RECURSIVE.stan"
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
# 21. STAN DATA
# =============================================================================

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


# =============================================================================
# 22. FIT FROZEN M4
# =============================================================================

cat(
  "\nFITTING FROZEN M4 AT WEEK 40\n"
)


fit <-
  model$sample(
    
    data =
      stan_data,
    
    seed =
      SEED,
    
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


saveRDS(
  fit,
  file.path(
    OUTPUT_DIR,
    "PACE_FLU_M4_WEEK40_FIT.rds"
  )
)


# =============================================================================
# 23. CONVERGENCE
# =============================================================================

conv_variables <- c(
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


sm <-
  fit$summary(
    variables =
      conv_variables
  )


diagnostic <-
  fit$diagnostic_summary()


convergence_summary <-
  tibble::tibble(
    
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
        diagnostic$num_divergent,
        na.rm = TRUE
      ),
    
    treedepth_hits =
      sum(
        diagnostic$num_max_treedepth,
        na.rm = TRUE
      )
  ) %>%
  
  dplyr::mutate(
    
    convergence_ok =
      
      .data$max_rhat <
      1.01 &
      
      .data$min_bulk_ESS >
      400 &
      
      .data$min_tail_ESS >
      400 &
      
      .data$divergences ==
      0 &
      
      .data$treedepth_hits ==
      0
  )


readr::write_csv(
  
  convergence_summary,
  
  file.path(
    OUTPUT_DIR,
    "01_MODEL_DIAGNOSTICS.csv"
  )
)


cat(
  "\nMODEL DIAGNOSTICS\n"
)

safe_print(
  convergence_summary
)


# =============================================================================
# 24. EXTRACT POSTERIOR PARAMETERS
# =============================================================================

parameter_variables <- c(
  "alpha",
  "alpha_state",
  "beta_memory",
  "beta_national",
  "beta_borrow",
  "beta_phase",
  "beta_interaction",
  "log_phi"
)


posterior_matrix <-
  fit$draws(
    
    variables =
      parameter_variables,
    
    format =
      "matrix"
  )


TOTAL_DRAWS <-
  nrow(
    posterior_matrix
  )


if (
  N_RECURSIVE_DRAWS >
  TOTAL_DRAWS
) {
  
  N_RECURSIVE_DRAWS <-
    TOTAL_DRAWS
}


draw_indices <-
  unique(
    round(
      seq(
        1,
        TOTAL_DRAWS,
        length.out =
          N_RECURSIVE_DRAWS
      )
    )
  )


posterior_matrix <-
  posterior_matrix[
    draw_indices,
    ,
    drop = FALSE
  ]


N_DRAWS <-
  nrow(
    posterior_matrix
  )


cat(
  "\nRecursive posterior trajectories: ",
  N_DRAWS,
  "\n",
  sep = ""
)


# =============================================================================
# 25. PARAMETER ARRAYS
# =============================================================================

alpha_draw <-
  posterior_matrix[
    ,
    "alpha"
  ]


beta_memory_draw <-
  posterior_matrix[
    ,
    "beta_memory"
  ]


beta_national_draw <-
  posterior_matrix[
    ,
    "beta_national"
  ]


beta_borrow_draw <-
  posterior_matrix[
    ,
    "beta_borrow"
  ]


beta_phase_draw <-
  posterior_matrix[
    ,
    "beta_phase"
  ]


beta_interaction_draw <-
  posterior_matrix[
    ,
    "beta_interaction"
  ]


phi_draw <-
  0.05 +
  exp(
    posterior_matrix[
      ,
      "log_phi"
    ]
  )


alpha_state_draw <-
  matrix(
    NA_real_,
    nrow = N_DRAWS,
    ncol = N_STATES
  )


for (
  state_idx in
  seq_len(N_STATES)
) {
  
  variable_name <-
    paste0(
      "alpha_state[",
      state_idx,
      "]"
    )
  
  alpha_state_draw[
    ,
    state_idx
  ] <-
    posterior_matrix[
      ,
      variable_name
    ]
}


# =============================================================================
# 26. OBSERVED HISTORY THROUGH WEEK 40
#
# Dimensions:
#   draw x state x week
#
# Weeks 1--40 are replicated observed history.
# Weeks 41--52 will be recursively simulated.
# =============================================================================

case_array <-
  array(
    
    NA_real_,
    
    dim = c(
      N_DRAWS,
      N_STATES,
      52L
    )
  )


observed_target <-
  dat %>%
  
  dplyr::filter(
    .data$season ==
      TARGET_SEASON
  ) %>%
  
  dplyr::arrange(
    .data$state_id,
    .data$week
  )


for (
  state_idx in
  seq_len(N_STATES)
) {
  
  state_history <-
    observed_target %>%
    
    dplyr::filter(
      .data$state_id ==
        state_idx
    ) %>%
    
    dplyr::arrange(
      .data$week
    )
  
  
  observed_cases_1_40 <-
    state_history$cases[
      state_history$week <=
        FORECAST_ORIGIN
    ]
  
  
  if (
    length(
      observed_cases_1_40
    ) !=
    FORECAST_ORIGIN
  ) {
    
    stop(
      paste0(
        "Incomplete observed history for state ",
        state_idx
      )
    )
  }
  
  
  case_array[
    ,
    state_idx,
    1:FORECAST_ORIGIN
  ] <-
    
    matrix(
      
      observed_cases_1_40,
      
      nrow =
        N_DRAWS,
      
      ncol =
        FORECAST_ORIGIN,
      
      byrow =
        TRUE
    )
}


# =============================================================================
# 27. POPULATION VECTOR
# =============================================================================

population_vector <-
  state_lookup$population


log_offset_vector <-
  log(
    population_vector /
      100000
  )


# =============================================================================
# 28. RECURSIVE SIMULATION
# =============================================================================

cat(
  "\nBEGINNING RECURSIVE FORECASTING...\n"
)


for (
  future_week in
  FORECAST_WEEKS
) {
  
  cat(
    "Simulating week ",
    future_week,
    " of 52\n",
    sep = ""
  )
  
  
  previous_week <-
    future_week - 1L
  
  
  for (
    draw_idx in
    seq_len(N_DRAWS)
  ) {
    
    # =========================================================================
    # Convert trajectory history to log incidence
    # =========================================================================
    
    log_incidence_history <-
      matrix(
        
        NA_real_,
        
        nrow =
          N_STATES,
        
        ncol =
          previous_week
      )
    
    
    for (
      state_idx in
      seq_len(N_STATES)
    ) {
      
      incidence_history <-
        100000 *
        case_array[
          draw_idx,
          state_idx,
          1:previous_week
        ] /
        population_vector[
          state_idx
        ]
      
      
      log_incidence_history[
        state_idx,
      ] <-
        log1p(
          incidence_history
        )
    }
    
    
    # =========================================================================
    # CURRENT LOG INCIDENCE AT t
    # =========================================================================
    
    current_log_incidence <-
      log_incidence_history[
        ,
        previous_week
      ]
    
    
    # =========================================================================
    # MEMORY AT t
    # =========================================================================
    
    if (
      previous_week <
      L_MEMORY
    ) {
      
      stop(
        "Insufficient memory history."
      )
    }
    
    
    memory_current <-
      rep(
        NA_real_,
        N_STATES
      )
    
    
    for (
      state_idx in
      seq_len(N_STATES)
    ) {
      
      recent_values <-
        log_incidence_history[
          state_idx,
          previous_week -
            (0:(L_MEMORY - 1L))
        ]
      
      
      memory_current[
        state_idx
      ] <-
        sum(
          MEMORY_WEIGHTS *
            recent_values
        )
    }
    
    
    # =========================================================================
    # PHASE AT t
    # =========================================================================
    
    previous_log_incidence <-
      log_incidence_history[
        ,
        previous_week - 1L
      ]
    
    
    growth_current <-
      current_log_incidence -
      previous_log_incidence
    
    
    phase_current <-
      plogis(
        growth_current /
          growth_scale
      )
    
    
    # =========================================================================
    # BORROWING AT t
    # =========================================================================
    
    borrowing_current <-
      as.numeric(
        W %*%
          current_log_incidence
      )
    
    
    borrow_phase_current <-
      borrowing_current *
      phase_current
    
    
    # =========================================================================
    # NATIONAL SIGNAL AT t
    # =========================================================================
    
    current_cases <-
      case_array[
        draw_idx,
        ,
        previous_week
      ]
    
    
    national_cases_current <-
      sum(
        current_cases
      )
    
    
    national_population <-
      sum(
        population_vector
      )
    
    
    national_incidence_current <-
      100000 *
      national_cases_current /
      national_population
    
    
    national_current <-
      log1p(
        national_incidence_current
      )
    
    
    # =========================================================================
    # LINEAR PREDICTOR FOR WEEK t+1
    # =========================================================================
    
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
    
    
    check_finite(
      eta,
      paste0(
        "recursive eta week ",
        future_week,
        " draw ",
        draw_idx
      )
    )
    
    
    # =========================================================================
    # DRAW NEXT-WEEK CASES
    #
    # R's rnbinom uses:
    #   mean = mu
    #   size = phi
    #
    # matching NB2:
    #   Var(Y) = mu + mu^2 / phi
    # =========================================================================
    
    mu <-
      exp(eta)
    
    
    if (
      any(
        !is.finite(mu)
      )
    ) {
      
      stop(
        paste0(
          "Non-finite mu at week ",
          future_week,
          ", draw ",
          draw_idx
        )
      )
    }
    
    
    simulated_cases <-
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
    
    
    case_array[
      draw_idx,
      ,
      future_week
    ] <-
      simulated_cases
  }
}


cat(
  "\nRecursive simulation complete.\n"
)


# =============================================================================
# 29. SAVE POSTERIOR TRAJECTORIES
# =============================================================================

saveRDS(
  
  list(
    
    case_array =
      case_array,
    
    forecast_weeks =
      FORECAST_WEEKS,
    
    state_lookup =
      state_lookup,
    
    W =
      W,
    
    growth_scale =
      growth_scale,
    
    posterior_draw_indices =
      draw_indices
  ),
  
  file.path(
    OUTPUT_DIR,
    "posterior_recursive_draws.rds"
  )
)


# =============================================================================
# 30. HOLDOUT TRUTH
#
# IMPORTANT:
# Only now do weeks 41--52 enter evaluation.
# =============================================================================

truth_state <-
  dat %>%
  
  dplyr::filter(
    
    .data$season ==
      TARGET_SEASON,
    
    .data$week %in%
      FORECAST_WEEKS
  ) %>%
  
  dplyr::select(
    .data$state_id,
    .data$state,
    .data$zone,
    .data$macro_region,
    .data$population,
    .data$week,
    .data$cases
  ) %>%
  
  dplyr::arrange(
    .data$state_id,
    .data$week
  )


# =============================================================================
# 31. STATE WEEKLY FORECASTS
# =============================================================================

state_weekly_results <-
  list()

counter <- 1L


for (
  week_value in
  FORECAST_WEEKS
) {
  
  for (
    state_idx in
    seq_len(N_STATES)
  ) {
    
    draws <-
      case_array[
        ,
        state_idx,
        week_value
      ]
    
    
    truth_row <-
      truth_state %>%
      
      dplyr::filter(
        
        .data$state_id ==
          state_idx,
        
        .data$week ==
          week_value
      )
    
    
    summary_i <-
      summarise_draws(
        
        draws,
        
        actual =
          truth_row$cases[
            1
          ]
      )
    
    
    state_weekly_results[[
      counter
    ]] <-
      
      tibble::tibble(
        
        season =
          TARGET_SEASON,
        
        forecast_origin =
          FORECAST_ORIGIN,
        
        week =
          week_value,
        
        horizon =
          week_value -
          FORECAST_ORIGIN,
        
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
          ],
        
        population =
          population_vector[
            state_idx
          ]
      ) %>%
      
      dplyr::bind_cols(
        summary_i
      )
    
    
    counter <-
      counter + 1L
  }
}


state_weekly_forecasts <-
  dplyr::bind_rows(
    state_weekly_results
  )


readr::write_csv(
  
  state_weekly_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "02_STATE_WEEKLY_FORECASTS.csv"
  )
)


# =============================================================================
# 32. FUNCTION FOR DRAW-BY-DRAW SPATIAL AGGREGATION
# =============================================================================

aggregate_weekly_draws <- function(
    grouping_vector,
    grouping_name
) {
  
  groups <-
    unique(
      grouping_vector
    )
  
  
  output <- list()
  counter <- 1L
  
  
  for (
    group_value in
    groups
  ) {
    
    state_indices <-
      which(
        grouping_vector ==
          group_value
      )
    
    
    for (
      week_value in
      FORECAST_WEEKS
    ) {
      
      draws <-
        rowSums(
          
          case_array[
            ,
            state_indices,
            week_value,
            drop = FALSE
          ][
            ,
            ,
            1
          ]
        )
      
      
      actual <-
        truth_state %>%
        
        dplyr::filter(
          
          .data$state_id %in%
            state_indices,
          
          .data$week ==
            week_value
        ) %>%
        
        dplyr::summarise(
          
          actual =
            sum(
              .data$cases
            )
        ) %>%
        
        dplyr::pull(
          .data$actual
        )
      
      
      summary_i <-
        summarise_draws(
          draws,
          actual
        )
      
      
      output[[
        counter
      ]] <-
        
        tibble::tibble(
          
          season =
            TARGET_SEASON,
          
          forecast_origin =
            FORECAST_ORIGIN,
          
          week =
            week_value,
          
          horizon =
            week_value -
            FORECAST_ORIGIN,
          
          spatial_level =
            grouping_name,
          
          location =
            group_value
        ) %>%
        
        dplyr::bind_cols(
          summary_i
        )
      
      
      counter <-
        counter + 1L
    }
  }
  
  
  dplyr::bind_rows(
    output
  )
}


# =============================================================================
# 33. ZONE FORECASTS
# =============================================================================

zone_weekly_forecasts <-
  aggregate_weekly_draws(
    
    grouping_vector =
      state_lookup$zone,
    
    grouping_name =
      "zone"
  )


readr::write_csv(
  
  zone_weekly_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "03_ZONE_WEEKLY_FORECASTS.csv"
  )
)


# =============================================================================
# 34. NORTH/SOUTH FORECASTS
# =============================================================================

north_south_weekly_forecasts <-
  aggregate_weekly_draws(
    
    grouping_vector =
      state_lookup$macro_region,
    
    grouping_name =
      "macro_region"
  )


readr::write_csv(
  
  north_south_weekly_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "04_NORTH_SOUTH_WEEKLY_FORECASTS.csv"
  )
)


# =============================================================================
# 35. NATIONAL WEEKLY FORECASTS
# =============================================================================

national_weekly_results <-
  list()


for (
  week_idx in
  seq_along(
    FORECAST_WEEKS
  )
) {
  
  week_value <-
    FORECAST_WEEKS[
      week_idx
    ]
  
  
  draws <-
    rowSums(
      case_array[
        ,
        ,
        week_value,
        drop = FALSE
      ][
        ,
        ,
        1
      ]
    )
  
  
  actual <-
    truth_state %>%
    
    dplyr::filter(
      .data$week ==
        week_value
    ) %>%
    
    dplyr::summarise(
      
      actual =
        sum(
          .data$cases
        )
    ) %>%
    
    dplyr::pull(
      .data$actual
    )
  
  
  national_weekly_results[[
    week_idx
  ]] <-
    
    tibble::tibble(
      
      season =
        TARGET_SEASON,
      
      forecast_origin =
        FORECAST_ORIGIN,
      
      week =
        week_value,
      
      horizon =
        week_value -
        FORECAST_ORIGIN,
      
      spatial_level =
        "national",
      
      location =
        "Nigeria"
    ) %>%
    
    dplyr::bind_cols(
      
      summarise_draws(
        draws,
        actual
      )
    )
}


national_weekly_forecasts <-
  dplyr::bind_rows(
    national_weekly_results
  )


readr::write_csv(
  
  national_weekly_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "05_NATIONAL_WEEKLY_FORECASTS.csv"
  )
)


# =============================================================================
# 36. WEEKLY VALIDATION BY HORIZON
# =============================================================================

weekly_validation_by_horizon <-
  state_weekly_forecasts %>%
  
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
    
    .groups = "drop"
  )


readr::write_csv(
  
  weekly_validation_by_horizon,
  
  file.path(
    OUTPUT_DIR,
    "06_STATE_WEEKLY_VALIDATION.csv"
  )
)


# =============================================================================
# 37. SPATIAL-LEVEL VALIDATION
# =============================================================================

state_spatial_summary <-
  state_weekly_forecasts %>%
  
  dplyr::summarise(
    
    spatial_level =
      "state",
    
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
      )
  )


zone_spatial_summary <-
  zone_weekly_forecasts %>%
  
  dplyr::summarise(
    
    spatial_level =
      "zone",
    
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
      )
  )


macro_spatial_summary <-
  north_south_weekly_forecasts %>%
  
  dplyr::summarise(
    
    spatial_level =
      "North/South",
    
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
      )
  )


national_spatial_summary <-
  national_weekly_forecasts %>%
  
  dplyr::summarise(
    
    spatial_level =
      "national",
    
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
      )
  )


spatial_validation <-
  dplyr::bind_rows(
    
    state_spatial_summary,
    
    zone_spatial_summary,
    
    macro_spatial_summary,
    
    national_spatial_summary
  )


readr::write_csv(
  
  spatial_validation,
  
  file.path(
    OUTPUT_DIR,
    "07_SPATIAL_LEVEL_VALIDATION.csv"
  )
)


# =============================================================================
# 38. STATE SEASONAL TARGETS
#
# Full season:
# observed weeks 1--40 +
# simulated weeks 41--52.
# =============================================================================

state_seasonal_results <-
  list()


for (
  state_idx in
  seq_len(N_STATES)
) {
  
  state_trajectories <-
    case_array[
      ,
      state_idx,
      ,
      drop = FALSE
    ][
      ,
      1,
    ]
  
  
  cumulative_draws <-
    rowSums(
      state_trajectories
    )
  
  
  peak_week_draws <-
    apply(
      
      state_trajectories,
      
      1,
      
      which.max
    )
  
  
  peak_cases_draws <-
    apply(
      
      state_trajectories,
      
      1,
      
      max
    )
  
  
  peak_incidence_draws <-
    100000 *
    peak_cases_draws /
    population_vector[
      state_idx
    ]
  
  
  actual_state <-
    observed_target %>%
    
    dplyr::filter(
      .data$state_id ==
        state_idx
    ) %>%
    
    dplyr::arrange(
      .data$week
    )
  
  
  actual_cumulative <-
    sum(
      actual_state$cases
    )
  
  
  actual_peak_index <-
    which.max(
      actual_state$cases
    )
  
  
  actual_peak_week <-
    actual_state$week[
      actual_peak_index
    ]
  
  
  actual_peak_cases <-
    max(
      actual_state$cases
    )
  
  
  actual_peak_incidence <-
    100000 *
    actual_peak_cases /
    population_vector[
      state_idx
    ]
  
  
  cumulative_summary <-
    summarise_draws(
      cumulative_draws,
      actual_cumulative
    )
  
  
  peak_incidence_summary <-
    summarise_draws(
      peak_incidence_draws,
      actual_peak_incidence
    )
  
  
  peak_week_quantiles <-
    stats::quantile(
      
      peak_week_draws,
      
      probs =
        c(
          0.05,
          0.25,
          0.50,
          0.75,
          0.95
        ),
      
      names = FALSE,
      
      type = 1
    )
  
  
  state_seasonal_results[[
    state_idx
  ]] <-
    
    tibble::tibble(
      
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
        ],
      
      actual_cumulative_cases =
        actual_cumulative,
      
      cumulative_q05 =
        cumulative_summary$q05,
      
      cumulative_q25 =
        cumulative_summary$q25,
      
      cumulative_q50 =
        cumulative_summary$q50,
      
      cumulative_q75 =
        cumulative_summary$q75,
      
      cumulative_q95 =
        cumulative_summary$q95,
      
      actual_peak_week =
        actual_peak_week,
      
      peak_week_q05 =
        peak_week_quantiles[1],
      
      peak_week_q25 =
        peak_week_quantiles[2],
      
      peak_week_q50 =
        peak_week_quantiles[3],
      
      peak_week_q75 =
        peak_week_quantiles[4],
      
      peak_week_q95 =
        peak_week_quantiles[5],
      
      probability_exact_peak_week =
        mean(
          peak_week_draws ==
            actual_peak_week
        ),
      
      probability_peak_within_1_week =
        mean(
          abs(
            peak_week_draws -
              actual_peak_week
          ) <=
            1
        ),
      
      actual_peak_incidence =
        actual_peak_incidence,
      
      peak_incidence_q05 =
        peak_incidence_summary$q05,
      
      peak_incidence_q25 =
        peak_incidence_summary$q25,
      
      peak_incidence_q50 =
        peak_incidence_summary$q50,
      
      peak_incidence_q75 =
        peak_incidence_summary$q75,
      
      peak_incidence_q95 =
        peak_incidence_summary$q95
    )
}


state_seasonal_targets <-
  dplyr::bind_rows(
    state_seasonal_results
  )


readr::write_csv(
  
  state_seasonal_targets,
  
  file.path(
    OUTPUT_DIR,
    "08_STATE_SEASONAL_TARGETS.csv"
  )
)


# =============================================================================
# 39. FUNCTION FOR AGGREGATED SEASONAL TARGETS
# =============================================================================

aggregate_seasonal_targets <- function(
    grouping_vector,
    spatial_level
) {
  
  groups <-
    unique(
      grouping_vector
    )
  
  
  output <- list()
  
  
  for (
    group_idx in
    seq_along(groups)
  ) {
    
    group_value <-
      groups[
        group_idx
      ]
    
    
    state_indices <-
      which(
        grouping_vector ==
          group_value
      )
    
    
    group_population <-
      sum(
        population_vector[
          state_indices
        ]
      )
    
    
    # Draw x week
    group_trajectory <-
      matrix(
        
        0,
        
        nrow =
          N_DRAWS,
        
        ncol =
          52L
      )
    
    
    for (
      state_idx in
      state_indices
    ) {
      
      group_trajectory <-
        group_trajectory +
        
        case_array[
          ,
          state_idx,
          ,
          drop = FALSE
        ][
          ,
          1,
        ]
    }
    
    
    cumulative_draws <-
      rowSums(
        group_trajectory
      )
    
    
    peak_week_draws <-
      apply(
        
        group_trajectory,
        
        1,
        
        which.max
      )
    
    
    peak_cases_draws <-
      apply(
        
        group_trajectory,
        
        1,
        
        max
      )
    
    
    peak_incidence_draws <-
      100000 *
      peak_cases_draws /
      group_population
    
    
    truth_group <-
      observed_target %>%
      
      dplyr::filter(
        .data$state_id %in%
          state_indices
      ) %>%
      
      dplyr::group_by(
        .data$week
      ) %>%
      
      dplyr::summarise(
        
        cases =
          sum(
            .data$cases
          ),
        
        .groups =
          "drop"
      ) %>%
      
      dplyr::arrange(
        .data$week
      )
    
    
    actual_cumulative <-
      sum(
        truth_group$cases
      )
    
    
    actual_peak_index <-
      which.max(
        truth_group$cases
      )
    
    
    actual_peak_week <-
      truth_group$week[
        actual_peak_index
      ]
    
    
    actual_peak_cases <-
      max(
        truth_group$cases
      )
    
    
    actual_peak_incidence <-
      100000 *
      actual_peak_cases /
      group_population
    
    
    cumulative_q <-
      stats::quantile(
        
        cumulative_draws,
        
        probs =
          c(
            0.05,
            0.25,
            0.50,
            0.75,
            0.95
          ),
        
        names = FALSE,
        
        type = 8
      )
    
    
    peak_week_q <-
      stats::quantile(
        
        peak_week_draws,
        
        probs =
          c(
            0.05,
            0.25,
            0.50,
            0.75,
            0.95
          ),
        
        names = FALSE,
        
        type = 1
      )
    
    
    peak_incidence_q <-
      stats::quantile(
        
        peak_incidence_draws,
        
        probs =
          c(
            0.05,
            0.25,
            0.50,
            0.75,
            0.95
          ),
        
        names = FALSE,
        
        type = 8
      )
    
    
    output[[
      group_idx
    ]] <-
      
      tibble::tibble(
        
        spatial_level =
          spatial_level,
        
        location =
          group_value,
        
        actual_cumulative_cases =
          actual_cumulative,
        
        cumulative_q05 =
          cumulative_q[1],
        
        cumulative_q25 =
          cumulative_q[2],
        
        cumulative_q50 =
          cumulative_q[3],
        
        cumulative_q75 =
          cumulative_q[4],
        
        cumulative_q95 =
          cumulative_q[5],
        
        actual_peak_week =
          actual_peak_week,
        
        peak_week_q05 =
          peak_week_q[1],
        
        peak_week_q25 =
          peak_week_q[2],
        
        peak_week_q50 =
          peak_week_q[3],
        
        peak_week_q75 =
          peak_week_q[4],
        
        peak_week_q95 =
          peak_week_q[5],
        
        probability_exact_peak_week =
          mean(
            peak_week_draws ==
              actual_peak_week
          ),
        
        probability_peak_within_1_week =
          mean(
            abs(
              peak_week_draws -
                actual_peak_week
            ) <=
              1
          ),
        
        actual_peak_incidence =
          actual_peak_incidence,
        
        peak_incidence_q05 =
          peak_incidence_q[1],
        
        peak_incidence_q25 =
          peak_incidence_q[2],
        
        peak_incidence_q50 =
          peak_incidence_q[3],
        
        peak_incidence_q75 =
          peak_incidence_q[4],
        
        peak_incidence_q95 =
          peak_incidence_q[5]
      )
  }
  
  
  dplyr::bind_rows(
    output
  )
}


# =============================================================================
# 40. ZONE SEASONAL TARGETS
# =============================================================================

zone_seasonal_targets <-
  aggregate_seasonal_targets(
    
    state_lookup$zone,
    
    "zone"
  )


readr::write_csv(
  
  zone_seasonal_targets,
  
  file.path(
    OUTPUT_DIR,
    "09_ZONE_SEASONAL_TARGETS.csv"
  )
)


# =============================================================================
# 41. NORTH/SOUTH SEASONAL TARGETS
# =============================================================================

north_south_seasonal_targets <-
  aggregate_seasonal_targets(
    
    state_lookup$macro_region,
    
    "macro_region"
  )


readr::write_csv(
  
  north_south_seasonal_targets,
  
  file.path(
    OUTPUT_DIR,
    "10_NORTH_SOUTH_SEASONAL_TARGETS.csv"
  )
)


# =============================================================================
# 42. NATIONAL SEASONAL TARGETS
# =============================================================================

national_grouping <-
  rep(
    "Nigeria",
    N_STATES
  )


national_seasonal_targets <-
  aggregate_seasonal_targets(
    
    national_grouping,
    
    "national"
  )


readr::write_csv(
  
  national_seasonal_targets,
  
  file.path(
    OUTPUT_DIR,
    "11_NATIONAL_SEASONAL_TARGETS.csv"
  )
)


# =============================================================================
# 43. NATIONAL PEAK-WEEK DISTRIBUTION
# =============================================================================

national_trajectory <-
  apply(
    case_array,
    c(1, 3),
    sum
  )


national_peak_week_draws <-
  apply(
    national_trajectory,
    1,
    which.max
  )


peak_week_distribution <-
  tibble::tibble(
    
    week =
      1:52,
    
    probability =
      vapply(
        
        1:52,
        
        function(w) {
          
          mean(
            national_peak_week_draws ==
              w
          )
        },
        
        numeric(1)
      )
  )


readr::write_csv(
  
  peak_week_distribution,
  
  file.path(
    OUTPUT_DIR,
    "12_PEAK_WEEK_DISTRIBUTIONS.csv"
  )
)


# =============================================================================
# 44. CALIBRATION
# =============================================================================

calibration <-
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
          state_weekly_forecasts$covered50
        ),
        
        mean(
          state_weekly_forecasts$covered80
        ),
        
        mean(
          state_weekly_forecasts$covered90
        )
      )
  ) %>%
  
  dplyr::mutate(
    
    deviation =
      .data$observed -
      .data$nominal
  )


readr::write_csv(
  
  calibration,
  
  file.path(
    OUTPUT_DIR,
    "13_CALIBRATION.csv"
  )
)


# =============================================================================
# 45. OVERALL RECURSIVE PERFORMANCE
# =============================================================================

overall_recursive_performance <-
  state_weekly_forecasts %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    horizons =
      dplyr::n_distinct(
        .data$horizon
      ),
    
    MAE =
      mean(
        .data$absolute_error
      ),
    
    median_absolute_error =
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
    
    median_CRPS =
      stats::median(
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
      )
  )


readr::write_csv(
  
  overall_recursive_performance,
  
  file.path(
    OUTPUT_DIR,
    "14_FINAL_RECURSIVE_PERFORMANCE.csv"
  )
)


# =============================================================================
# 46. NATIONAL FORECAST FIGURE
# =============================================================================

national_observed <-
  observed_target %>%
  
  dplyr::group_by(
    .data$week
  ) %>%
  
  dplyr::summarise(
    
    cases =
      sum(
        .data$cases
      ),
    
    .groups = "drop"
  )


national_plot_data <-
  national_weekly_forecasts %>%
  
  dplyr::select(
    .data$week,
    .data$q05,
    .data$q25,
    .data$q50,
    .data$q75,
    .data$q95
  )


p_national <-
  ggplot2::ggplot() +
  
  ggplot2::geom_line(
    
    data =
      national_observed,
    
    ggplot2::aes(
      x = week,
      y = cases
    )
  ) +
  
  ggplot2::geom_ribbon(
    
    data =
      national_plot_data,
    
    ggplot2::aes(
      x = week,
      ymin = q05,
      ymax = q95
    ),
    
    alpha = 0.20
  ) +
  
  ggplot2::geom_ribbon(
    
    data =
      national_plot_data,
    
    ggplot2::aes(
      x = week,
      ymin = q25,
      ymax = q75
    ),
    
    alpha = 0.30
  ) +
  
  ggplot2::geom_line(
    
    data =
      national_plot_data,
    
    ggplot2::aes(
      x = week,
      y = q50
    ),
    
    linetype = 2
  ) +
  
  ggplot2::geom_vline(
    
    xintercept =
      FORECAST_ORIGIN,
    
    linetype = 3
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU Recursive National Influenza Forecast",
    
    subtitle =
      "Observed through week 40; recursive forecast weeks 41-52",
    
    x =
      "Epidemiological week",
    
    y =
      "Weekly cases"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "01_NATIONAL_RECURSIVE_FORECAST.png"
  ),
  
  p_national,
  
  width = 9,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 47. HORIZON PERFORMANCE FIGURE
# =============================================================================

p_horizon <-
  ggplot2::ggplot(
    
    weekly_validation_by_horizon,
    
    ggplot2::aes(
      x = horizon,
      y = mean_WIS
    )
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU Probabilistic Accuracy by Forecast Horizon",
    
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
  
  p_horizon,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 48. NATIONAL PEAK DISTRIBUTION FIGURE
# =============================================================================

p_peak <-
  ggplot2::ggplot(
    
    peak_week_distribution,
    
    ggplot2::aes(
      x = week,
      y = probability
    )
  ) +
  
  ggplot2::geom_col() +
  
  ggplot2::labs(
    
    title =
      "Posterior Distribution of National Influenza Peak Week",
    
    x =
      "Epidemiological week",
    
    y =
      "Posterior probability"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "03_NATIONAL_PEAK_WEEK_DISTRIBUTION.png"
  ),
  
  p_peak,
  
  width = 9,
  
  height = 6,
  
  dpi = 300
)


# =============================================================================
# 49. SESSION INFO
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
# 50. FINAL REPORT
# =============================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "PACE-FLU RECURSIVE CASE FORECAST COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nMODEL DIAGNOSTICS\n"
)

safe_print(
  convergence_summary
)


cat(
  "\nOVERALL 12-WEEK STATE-LEVEL PERFORMANCE\n"
)

safe_print(
  overall_recursive_performance
)


cat(
  "\nPERFORMANCE BY HORIZON\n"
)

safe_print(
  weekly_validation_by_horizon
)


cat(
  "\nSPATIAL PERFORMANCE\n"
)

safe_print(
  spatial_validation
)


cat(
  "\nCALIBRATION\n"
)

safe_print(
  calibration
)


cat(
  "\nNATIONAL SEASONAL TARGETS\n"
)

safe_print(
  national_seasonal_targets
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
  "END OF FINAL RECURSIVE CASE PIPELINE\n"
)

cat(
  "============================================================\n"
)