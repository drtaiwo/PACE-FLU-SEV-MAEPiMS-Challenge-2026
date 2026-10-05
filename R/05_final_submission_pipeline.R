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

# Final severity convergence safeguard (sampling only; model unchanged)
FINAL_RHAT_THRESHOLD <- 1.01
SEVERITY_REFIT_WARMUP <- 2000L
SEVERITY_REFIT_SAMPLING <- 3000L
SEVERITY_REFIT_ADAPT_DELTA <- 0.995


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
      ),
    
    fallback_targets =
      sum(vapply(seq_len(N_STATES), function(target_idx) {
        target_edges <- edges %>% dplyr::filter(.data$target_id == target_idx)
        weights <- ifelse(is.finite(target_edges$coefficient),
                          pmax(target_edges$coefficient, 0), 0)
        sum(weights) <= 0
      }, logical(1)))
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

# =============================================================================
# 47. FINAL INTEGRATED SEVERITY SUBMISSION EXTENSION
# =============================================================================
#
# This section extends the already-fitted frozen M4 recursive CASE forecast
# above to the validated severity cascade:
#
#   CASES -> HOSPITALIZATIONS -> DEATHS
#
# It preserves the case trajectories already generated by frozen M4.
# Hospitalization and mortality models use exactly the architecture validated
# in Script 04.  Uncertainty is propagated draw-by-draw:
#
#   C^(m)_{i,t+1} -> H^(m)_{i,t+1} -> D^(m)_{i,t+1}.
#
# IMPORTANT:
#   * Weeks 41--52 are not used for fitting.
#   * Previous seasons are fully available.
#   * In 2025/26, severity training predictor rows stop at week 39, so their
#     next-week outcomes stop at observed week 40.
#   * The case model is NOT changed here.
#   * Geographic aggregation is performed draw-by-draw; state quantiles are
#     never summed.
# =============================================================================

INTEGRATED_OUTPUT_DIR <- "PACE_FLU_SEV_FINAL_SUBMISSION"

dir.create(
  INTEGRATED_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

INTEGRATED_STAN_DIR <- file.path(INTEGRATED_OUTPUT_DIR, "stan")
dir.create(INTEGRATED_STAN_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 47.1 Add severity next-week outcomes without altering frozen case features
# -----------------------------------------------------------------------------

severity_dat <- dat %>%
  dplyr::group_by(.data$season, .data$state) %>%
  dplyr::arrange(.data$week, .by_group = TRUE) %>%
  dplyr::mutate(
    hospitalization_next = dplyr::lead(.data$hospitalizations, 1L),
    deaths_next = dplyr::lead(.data$deaths, 1L),
    log_hosp_current = log1p(.data$hospitalizations),
    log_death_current = log1p(.data$deaths)
  ) %>%
  dplyr::ungroup()

severity_train <- severity_dat %>%
  dplyr::filter(
    !is.na(.data$cases_next),
    !is.na(.data$hospitalization_next),
    !is.na(.data$deaths_next),
    !is.na(.data$growth_raw),
    (.data$season_index < TARGET_SEASON_INDEX) |
      (.data$season_index == TARGET_SEASON_INDEX &
         .data$week <= FORECAST_ORIGIN - 1L)
  ) %>%
  dplyr::mutate(
    p_growth = plogis(.data$growth_raw / growth_scale)
  ) %>%
  dplyr::filter(
    stats::complete.cases(
      .data$cases_next,
      .data$hospitalization_next,
      .data$deaths_next,
      .data$log_hosp_current,
      .data$log_death_current,
      .data$p_growth
    )
  )

if (nrow(severity_train) < 200L) {
  stop("Insufficient leakage-free severity training data.")
}

# -----------------------------------------------------------------------------
# 47.2 Validated hospitalization and mortality Stan models
# -----------------------------------------------------------------------------

hospital_stan_code <- "\ndata {\n\n  int<lower=1> N;\n\n  int<lower=1> J;\n\n  array[N]\n    int<lower=0>\n    y;\n\n  array[N]\n    int<lower=1,upper=J>\n    state_id;\n\n  vector[N]\n    log_case_next;\n\n  vector[N]\n    log_hosp_current;\n\n  vector[N]\n    phase;\n}\n\n\nparameters {\n\n  real alpha;\n\n  vector[J]\n    state_raw;\n\n  real<lower=0>\n    sigma_state;\n\n  real beta_case;\n\n  real beta_hosp_memory;\n\n  real beta_phase;\n\n  real log_phi;\n}\n\n\ntransformed parameters {\n\n  vector[J]\n    alpha_state;\n\n  alpha_state =\n    sigma_state *\n    state_raw;\n}\n\n\nmodel {\n\n  vector[N]\n    eta;\n\n  real phi;\n\n\n  alpha ~\n    normal(\n      -4,\n      2\n    );\n\n\n  state_raw ~\n    normal(\n      0,\n      1\n    );\n\n\n  sigma_state ~\n    normal(\n      0,\n      1\n    );\n\n\n  beta_case ~\n    normal(\n      1,\n      0.5\n    );\n\n\n  beta_hosp_memory ~\n    normal(\n      0,\n      1\n    );\n\n\n  beta_phase ~\n    normal(\n      0,\n      1\n    );\n\n\n  log_phi ~\n    normal(\n      log(5),\n      0.75\n    );\n\n\n  phi =\n    0.05 +\n    exp(\n      log_phi\n    );\n\n\n  eta =\n\n      alpha\n\n      + alpha_state[\n          state_id\n        ]\n\n      + beta_case *\n        log_case_next\n\n      + beta_hosp_memory *\n        log_hosp_current\n\n      + beta_phase *\n        phase;\n\n\n  y ~\n    neg_binomial_2_log(\n      eta,\n      phi\n    );\n}"

death_stan_code <- "\ndata {\n\n  int<lower=1> N;\n\n  int<lower=1> J;\n\n  array[N]\n    int<lower=0>\n    y;\n\n  array[N]\n    int<lower=1,upper=J>\n    state_id;\n\n  vector[N]\n    log_case_next;\n\n  vector[N]\n    log_hosp_next;\n\n  vector[N]\n    log_death_current;\n\n  vector[N]\n    phase;\n}\n\n\nparameters {\n\n  real alpha;\n\n  vector[J]\n    state_raw;\n\n  real<lower=0>\n    sigma_state;\n\n  real beta_case;\n\n  real beta_hosp;\n\n  real beta_death_memory;\n\n  real beta_phase;\n\n  real log_phi;\n}\n\n\ntransformed parameters {\n\n  vector[J]\n    alpha_state;\n\n  alpha_state =\n    sigma_state *\n    state_raw;\n}\n\n\nmodel {\n\n  vector[N]\n    eta;\n\n  real phi;\n\n\n  alpha ~\n    normal(\n      -5,\n      2\n    );\n\n\n  state_raw ~\n    normal(\n      0,\n      1\n    );\n\n\n  sigma_state ~\n    normal(\n      0,\n      1\n    );\n\n\n  beta_case ~\n    normal(\n      0.5,\n      0.5\n    );\n\n\n  beta_hosp ~\n    normal(\n      0.5,\n      0.5\n    );\n\n\n  beta_death_memory ~\n    normal(\n      0,\n      1\n    );\n\n\n  beta_phase ~\n    normal(\n      0,\n      1\n    );\n\n\n  log_phi ~\n    normal(\n      log(5),\n      0.75\n    );\n\n\n  phi =\n    0.05 +\n    exp(\n      log_phi\n    );\n\n\n  eta =\n\n      alpha\n\n      + alpha_state[\n          state_id\n        ]\n\n      + beta_case *\n        log_case_next\n\n      + beta_hosp *\n        log_hosp_next\n\n      + beta_death_memory *\n        log_death_current\n\n      + beta_phase *\n        phase;\n\n\n  y ~\n    neg_binomial_2_log(\n      eta,\n      phi\n    );\n}"

HOSPITAL_STAN_FILE <- file.path(
  INTEGRATED_STAN_DIR,
  "PACE_FLU_HOSPITALIZATION_MODEL.stan"
)

DEATH_STAN_FILE <- file.path(
  INTEGRATED_STAN_DIR,
  "PACE_FLU_MORTALITY_MODEL.stan"
)

writeLines(hospital_stan_code, HOSPITAL_STAN_FILE)
writeLines(death_stan_code, DEATH_STAN_FILE)

cat("\nCompiling final hospitalization model...\n")
hospital_model_final <- cmdstanr::cmdstan_model(HOSPITAL_STAN_FILE)

cat("\nCompiling final mortality model...\n")
death_model_final <- cmdstanr::cmdstan_model(DEATH_STAN_FILE)

# -----------------------------------------------------------------------------
# 47.3 Stan data
# -----------------------------------------------------------------------------

hospital_data_final <- list(
  N = nrow(severity_train),
  J = N_STATES,
  y = as.integer(severity_train$hospitalization_next),
  state_id = as.integer(severity_train$state_id),
  log_case_next = log1p(severity_train$cases_next),
  log_hosp_current = as.vector(severity_train$log_hosp_current),
  phase = as.vector(severity_train$p_growth)
)

death_data_final <- list(
  N = nrow(severity_train),
  J = N_STATES,
  y = as.integer(severity_train$deaths_next),
  state_id = as.integer(severity_train$state_id),
  log_case_next = log1p(severity_train$cases_next),
  log_hosp_next = log1p(severity_train$hospitalization_next),
  log_death_current = as.vector(severity_train$log_death_current),
  phase = as.vector(severity_train$p_growth)
)

# -----------------------------------------------------------------------------
# 47.4 Fit severity models
#
# Mortality gets stricter sampling because Script 04 identified a few death
# fits with Rhat around 1.02.  If max Rhat remains > 1.01, refit automatically.
# -----------------------------------------------------------------------------

cat("\nFitting final hospitalization model...\n")

fit_hospital_final <- hospital_model_final$sample(
  data = hospital_data_final,
  seed = SEED + 300001L,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  step_size = 0.01,
  init = 0,
  refresh = 500
)

cat("\nFitting final mortality model...\n")

fit_death_final <- death_model_final$sample(
  data = death_data_final,
  seed = SEED + 400001L,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 1500L,
  iter_sampling = 2000L,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  step_size = 0.01,
  init = 0,
  refresh = 500
)

get_final_diagnostics <- function(fit, model_name) {
  sm <- fit$summary()
  ds <- fit$diagnostic_summary()
  tibble::tibble(
    model = model_name,
    max_rhat = max(sm$rhat, na.rm = TRUE),
    min_bulk_ESS = min(sm$ess_bulk, na.rm = TRUE),
    min_tail_ESS = min(sm$ess_tail, na.rm = TRUE),
    divergences = sum(ds$num_divergent, na.rm = TRUE),
    treedepth_hits = sum(ds$num_max_treedepth, na.rm = TRUE)
  )
}

# Automatic convergence safeguards. These refits change sampling effort only.
hospital_diag_first <- get_final_diagnostics(fit_hospital_final, "Hospitalizations")
cat("\nInitial hospitalization diagnostics:\n")
print(hospital_diag_first, n = Inf)

if (is.finite(hospital_diag_first$max_rhat) &&
    hospital_diag_first$max_rhat > FINAL_RHAT_THRESHOLD) {
  cat("\nHospitalization max Rhat = ", round(hospital_diag_first$max_rhat, 6),
      " > ", FINAL_RHAT_THRESHOLD,
      ". Refitting with stricter sampling...\n", sep = "")
  fit_hospital_final <- hospital_model_final$sample(
    data = hospital_data_final, seed = SEED + 300002L,
    chains = 4L, parallel_chains = 4L,
    iter_warmup = SEVERITY_REFIT_WARMUP,
    iter_sampling = SEVERITY_REFIT_SAMPLING,
    adapt_delta = SEVERITY_REFIT_ADAPT_DELTA,
    max_treedepth = 15L, step_size = 0.01, init = 0, refresh = 500
  )
  cat("\nHospitalization diagnostics after refit:\n")
  print(get_final_diagnostics(fit_hospital_final, "Hospitalizations"), n = Inf)
} else {
  cat("\nHospitalization model passed max Rhat <= ", FINAL_RHAT_THRESHOLD,
      "; no refit required.\n", sep = "")
}

death_diag_first <- get_final_diagnostics(fit_death_final, "Deaths")
cat("\nInitial mortality diagnostics:\n")
print(death_diag_first, n = Inf)

if (is.finite(death_diag_first$max_rhat) &&
    death_diag_first$max_rhat > FINAL_RHAT_THRESHOLD) {
  cat("\nMortality max Rhat = ", round(death_diag_first$max_rhat, 6),
      " > ", FINAL_RHAT_THRESHOLD,
      ". Refitting with stricter sampling...\n", sep = "")
  fit_death_final <- death_model_final$sample(
    data = death_data_final, seed = SEED + 400002L,
    chains = 4L, parallel_chains = 4L,
    iter_warmup = SEVERITY_REFIT_WARMUP,
    iter_sampling = SEVERITY_REFIT_SAMPLING,
    adapt_delta = SEVERITY_REFIT_ADAPT_DELTA,
    max_treedepth = 15L, step_size = 0.01, init = 0, refresh = 500
  )
  cat("\nMortality diagnostics after refit:\n")
  print(get_final_diagnostics(fit_death_final, "Deaths"), n = Inf)
} else {
  cat("\nMortality model passed max Rhat <= ", FINAL_RHAT_THRESHOLD,
      "; no refit required.\n", sep = "")
}

final_diagnostics <- dplyr::bind_rows(
  get_final_diagnostics(fit, "Cases"),
  get_final_diagnostics(fit_hospital_final, "Hospitalizations"),
  get_final_diagnostics(fit_death_final, "Deaths")
) %>%
  dplyr::mutate(
    convergence_ok = .data$max_rhat <= FINAL_RHAT_THRESHOLD &
      .data$min_bulk_ESS > 400 & .data$min_tail_ESS > 400 &
      .data$divergences == 0 & .data$treedepth_hits == 0
  )

readr::write_csv(
  final_diagnostics,
  file.path(INTEGRATED_OUTPUT_DIR, "01_FINAL_MODEL_DIAGNOSTICS.csv")
)

if (any(!final_diagnostics$convergence_ok)) {
  warning("At least one final model failed the convergence criteria; inspect 01_FINAL_MODEL_DIAGNOSTICS.csv before submission.")
}

# -----------------------------------------------------------------------------
# 47.5 Extract severity posterior draws and align to case trajectories
# -----------------------------------------------------------------------------

extract_selected_draws <- function(fit_object, variables, n_keep) {
  
  mat <- fit_object$draws(
    variables = variables,
    format = "matrix"
  )
  
  idx <- unique(round(seq(1, nrow(mat), length.out = n_keep)))
  
  mat[idx, , drop = FALSE]
}

hospital_draws_final <- extract_selected_draws(
  fit_hospital_final,
  c(
    "alpha",
    "alpha_state",
    "beta_case",
    "beta_hosp_memory",
    "beta_phase",
    "log_phi"
  ),
  N_DRAWS
)

death_draws_final <- extract_selected_draws(
  fit_death_final,
  c(
    "alpha",
    "alpha_state",
    "beta_case",
    "beta_hosp",
    "beta_death_memory",
    "beta_phase",
    "log_phi"
  ),
  N_DRAWS
)

# Defensive alignment if posterior extraction produced a slightly different
# number of unique rows.
N_INTEGRATED_DRAWS <- min(
  N_DRAWS,
  nrow(hospital_draws_final),
  nrow(death_draws_final)
)

case_array_integrated <- case_array[
  seq_len(N_INTEGRATED_DRAWS),
  ,
  ,
  drop = FALSE
]

hospital_draws_final <- hospital_draws_final[
  seq_len(N_INTEGRATED_DRAWS),
  ,
  drop = FALSE
]

death_draws_final <- death_draws_final[
  seq_len(N_INTEGRATED_DRAWS),
  ,
  drop = FALSE
]

extract_state_effect_matrix <- function(draw_matrix, n_states) {
  
  out <- matrix(
    NA_real_,
    nrow = nrow(draw_matrix),
    ncol = n_states
  )
  
  for (j in seq_len(n_states)) {
    nm <- paste0("alpha_state[", j, "]")
    out[, j] <- draw_matrix[, nm]
  }
  
  out
}

hospital_alpha_state <- extract_state_effect_matrix(
  hospital_draws_final,
  N_STATES
)

death_alpha_state <- extract_state_effect_matrix(
  death_draws_final,
  N_STATES
)

hospital_phi <- 0.05 + exp(hospital_draws_final[, "log_phi"])
death_phi <- 0.05 + exp(death_draws_final[, "log_phi"])

# -----------------------------------------------------------------------------
# 47.6 Initialize observed hospitalization/death histories through week 40
# -----------------------------------------------------------------------------

hospital_array <- array(
  NA_real_,
  dim = c(N_INTEGRATED_DRAWS, N_STATES, 52L)
)

death_array <- array(
  NA_real_,
  dim = c(N_INTEGRATED_DRAWS, N_STATES, 52L)
)

for (state_idx in seq_len(N_STATES)) {
  
  hist_state <- severity_dat %>%
    dplyr::filter(
      .data$season == TARGET_SEASON,
      .data$state_id == state_idx,
      .data$week <= FORECAST_ORIGIN
    ) %>%
    dplyr::arrange(.data$week)
  
  if (nrow(hist_state) != FORECAST_ORIGIN) {
    stop(
      paste0(
        "Incomplete severity history through week 40 for state ",
        state_idx
      )
    )
  }
  
  hospital_array[, state_idx, 1:FORECAST_ORIGIN] <- matrix(
    hist_state$hospitalizations,
    nrow = N_INTEGRATED_DRAWS,
    ncol = FORECAST_ORIGIN,
    byrow = TRUE
  )
  
  death_array[, state_idx, 1:FORECAST_ORIGIN] <- matrix(
    hist_state$deaths,
    nrow = N_INTEGRATED_DRAWS,
    ncol = FORECAST_ORIGIN,
    byrow = TRUE
  )
}

# -----------------------------------------------------------------------------
# 47.7 Recursive severity propagation using frozen M4 case trajectories
# -----------------------------------------------------------------------------

cat("\nBEGINNING INTEGRATED C -> H -> D RECURSIVE PROPAGATION...\n")

for (future_week in FORECAST_WEEKS) {
  
  previous_week <- future_week - 1L
  
  cat(
    "Propagating severity for week ",
    future_week,
    " of 52\n",
    sep = ""
  )
  
  for (draw_idx in seq_len(N_INTEGRATED_DRAWS)) {
    
    # Epidemic phase is recomputed from this draw's recursive case trajectory.
    current_incidence <- 100000 *
      case_array_integrated[draw_idx, , previous_week] /
      population_vector
    
    previous_incidence <- 100000 *
      case_array_integrated[draw_idx, , previous_week - 1L] /
      population_vector
    
    current_log_incidence <- log1p(current_incidence)
    previous_log_incidence <- log1p(previous_incidence)
    
    growth_current <- current_log_incidence - previous_log_incidence
    phase_current <- plogis(growth_current / growth_scale)
    
    predicted_cases <- case_array_integrated[
      draw_idx,
      ,
      future_week
    ]
    
    current_hospitalizations <- hospital_array[
      draw_idx,
      ,
      previous_week
    ]
    
    eta_hospital <-
      hospital_draws_final[draw_idx, "alpha"] +
      hospital_alpha_state[draw_idx, ] +
      hospital_draws_final[draw_idx, "beta_case"] *
      log1p(predicted_cases) +
      hospital_draws_final[draw_idx, "beta_hosp_memory"] *
      log1p(current_hospitalizations) +
      hospital_draws_final[draw_idx, "beta_phase"] *
      phase_current
    
    mu_hospital <- exp(eta_hospital)
    
    if (any(!is.finite(mu_hospital))) {
      stop(
        paste0(
          "Non-finite hospitalization mean at week ",
          future_week,
          ", draw ",
          draw_idx
        )
      )
    }
    
    predicted_hospitalizations <- stats::rnbinom(
      N_STATES,
      mu = mu_hospital,
      size = hospital_phi[draw_idx]
    )
    
    hospital_array[
      draw_idx,
      ,
      future_week
    ] <- predicted_hospitalizations
    
    current_deaths <- death_array[
      draw_idx,
      ,
      previous_week
    ]
    
    eta_death <-
      death_draws_final[draw_idx, "alpha"] +
      death_alpha_state[draw_idx, ] +
      death_draws_final[draw_idx, "beta_case"] *
      log1p(predicted_cases) +
      death_draws_final[draw_idx, "beta_hosp"] *
      log1p(predicted_hospitalizations) +
      death_draws_final[draw_idx, "beta_death_memory"] *
      log1p(current_deaths) +
      death_draws_final[draw_idx, "beta_phase"] *
      phase_current
    
    mu_death <- exp(eta_death)
    
    if (any(!is.finite(mu_death))) {
      stop(
        paste0(
          "Non-finite death mean at week ",
          future_week,
          ", draw ",
          draw_idx
        )
      )
    }
    
    predicted_deaths <- stats::rnbinom(
      N_STATES,
      mu = mu_death,
      size = death_phi[draw_idx]
    )
    
    death_array[
      draw_idx,
      ,
      future_week
    ] <- predicted_deaths
  }
}

cat("\nIntegrated recursive propagation complete.\n")

# -----------------------------------------------------------------------------
# 47.8 Quantile helpers
# -----------------------------------------------------------------------------

FINAL_PROBS <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
FINAL_QNAMES <- c("q05", "q10", "q25", "median", "q75", "q90", "q95")

summarise_vector_quantiles <- function(x) {
  
  q <- stats::quantile(
    x,
    probs = FINAL_PROBS,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
  
  stats::setNames(as.list(as.numeric(q)), FINAL_QNAMES)
}

array_for_outcome <- list(
  Cases = case_array_integrated,
  Hospitalizations = hospital_array,
  Deaths = death_array
)

# -----------------------------------------------------------------------------
# 47.9 Draw-wise geographic aggregation
# -----------------------------------------------------------------------------

state_lookup_final <- state_lookup %>%
  dplyr::mutate(
    macro_region = dplyr::case_when(
      .data$zone %in% c("NW", "NE", "NC") ~ "North",
      .data$zone %in% c("SW", "SE", "SS") ~ "South",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(state_lookup_final$macro_region))) {
  stop("Could not assign North/South macro-region to all states.")
}

aggregate_indices <- list()

# States
for (j in seq_len(N_STATES)) {
  aggregate_indices[[paste0("State||", state_lookup_final$state[j])]] <- j
}

# Zones
for (z in sort(unique(state_lookup_final$zone))) {
  aggregate_indices[[paste0("Zone||", z)]] <-
    which(state_lookup_final$zone == z)
}

# North/South
for (m in c("North", "South")) {
  aggregate_indices[[paste0("Macro||", m)]] <-
    which(state_lookup_final$macro_region == m)
}

# National
aggregate_indices[["National||Nigeria"]] <- seq_len(N_STATES)

population_by_unit <- function(indices) {
  sum(population_vector[indices])
}

weekly_rows <- list()
weekly_counter <- 1L

for (outcome_name in names(array_for_outcome)) {
  
  arr <- array_for_outcome[[outcome_name]]
  
  for (future_week in FORECAST_WEEKS) {
    
    for (unit_key in names(aggregate_indices)) {
      
      parts <- strsplit(unit_key, "\\|\\|")[[1]]
      level_name <- parts[1]
      location_name <- parts[2]
      idx <- aggregate_indices[[unit_key]]
      
      # Draw-wise aggregation: sum states inside each posterior trajectory.
      draw_counts <- rowSums(
        arr[
          ,
          idx,
          future_week,
          drop = FALSE
        ][, , 1, drop = FALSE],
        na.rm = TRUE
      )
      
      unit_population <- population_by_unit(idx)
      
      draw_rate <- 100000 * draw_counts / unit_population
      
      q_count <- summarise_vector_quantiles(draw_counts)
      q_rate <- summarise_vector_quantiles(draw_rate)
      
      row <- tibble::tibble(
        model = "PACE-FLU-SEV",
        season = TARGET_SEASON,
        origin_week = FORECAST_ORIGIN,
        target_week = future_week,
        horizon = future_week - FORECAST_ORIGIN,
        outcome = outcome_name,
        spatial_level = level_name,
        location = location_name,
        population = unit_population
      )
      
      for (nm in FINAL_QNAMES) {
        row[[nm]] <- q_count[[nm]]
        row[[paste0(nm, "_per100k")]] <- q_rate[[nm]]
      }
      
      weekly_rows[[weekly_counter]] <- row
      weekly_counter <- weekly_counter + 1L
    }
  }
}

weekly_submission_wide <- dplyr::bind_rows(weekly_rows) %>%
  dplyr::arrange(
    factor(.data$outcome, levels = c("Cases", "Hospitalizations", "Deaths")),
    .data$target_week,
    .data$spatial_level,
    .data$location
  )

# Expected:
# 3 outcomes x 12 weeks x (37 states + 6 zones + 2 macro + 1 national)
EXPECTED_WEEKLY_ROWS <- 3L * 12L * (37L + 6L + 2L + 1L)

if (nrow(weekly_submission_wide) != EXPECTED_WEEKLY_ROWS) {
  stop(
    paste0(
      "Unexpected number of weekly submission rows: ",
      nrow(weekly_submission_wide),
      "; expected ",
      EXPECTED_WEEKLY_ROWS
    )
  )
}

# -----------------------------------------------------------------------------
# 47.10 Direction probabilities
# -----------------------------------------------------------------------------

direction_rows <- list()
direction_counter <- 1L

for (outcome_name in names(array_for_outcome)) {
  
  arr <- array_for_outcome[[outcome_name]]
  
  for (future_week in FORECAST_WEEKS) {
    
    previous_week <- future_week - 1L
    
    for (state_idx in seq_len(N_STATES)) {
      
      current_draw <- arr[, state_idx, future_week]
      previous_draw <- arr[, state_idx, previous_week]
      
      direction_rows[[direction_counter]] <- tibble::tibble(
        season = TARGET_SEASON,
        origin_week = FORECAST_ORIGIN,
        target_week = future_week,
        horizon = future_week - FORECAST_ORIGIN,
        outcome = outcome_name,
        state = state_lookup_final$state[state_idx],
        probability_increasing = mean(current_draw > previous_draw),
        probability_decreasing = mean(current_draw < previous_draw),
        probability_unchanged = mean(current_draw == previous_draw)
      )
      
      direction_counter <- direction_counter + 1L
    }
  }
}

direction_probabilities <- dplyr::bind_rows(direction_rows)

# -----------------------------------------------------------------------------
# 47.11 Seasonal cumulative targets and case attack rate
#
# These combine OBSERVED weeks 1--40 with simulated weeks 41--52.
# They are not purely prospective full-season forecasts.
# -----------------------------------------------------------------------------

seasonal_rows <- list()
seasonal_counter <- 1L

for (outcome_name in names(array_for_outcome)) {
  
  arr <- array_for_outcome[[outcome_name]]
  
  for (unit_key in names(aggregate_indices)) {
    
    parts <- strsplit(unit_key, "\\|\\|")[[1]]
    level_name <- parts[1]
    location_name <- parts[2]
    idx <- aggregate_indices[[unit_key]]
    
    cumulative_draws <- rep(0, N_INTEGRATED_DRAWS)
    
    for (wk in 1:52) {
      cumulative_draws <- cumulative_draws +
        rowSums(
          arr[, idx, wk, drop = FALSE][, , 1, drop = FALSE],
          na.rm = TRUE
        )
    }
    
    unit_population <- population_by_unit(idx)
    q_cum <- summarise_vector_quantiles(cumulative_draws)
    
    row <- tibble::tibble(
      model = "PACE-FLU-SEV",
      season = TARGET_SEASON,
      origin_week = FORECAST_ORIGIN,
      outcome = outcome_name,
      spatial_level = level_name,
      location = location_name,
      population = unit_population
    )
    
    for (nm in FINAL_QNAMES) {
      row[[paste0("cumulative_", nm)]] <- q_cum[[nm]]
    }
    
    if (outcome_name == "Cases") {
      
      attack_draws <- 100 * cumulative_draws / unit_population
      q_attack <- summarise_vector_quantiles(attack_draws)
      
      for (nm in FINAL_QNAMES) {
        row[[paste0("attack_rate_percent_", nm)]] <- q_attack[[nm]]
      }
      
    } else {
      
      for (nm in FINAL_QNAMES) {
        row[[paste0("attack_rate_percent_", nm)]] <- NA_real_
      }
    }
    
    seasonal_rows[[seasonal_counter]] <- row
    seasonal_counter <- seasonal_counter + 1L
  }
}

seasonal_targets <- dplyr::bind_rows(seasonal_rows)

# -----------------------------------------------------------------------------
# 47.12 Peak week and peak incidence distributions for CASES
# -----------------------------------------------------------------------------

peak_rows <- list()
peak_counter <- 1L
peak_summary_rows <- list()
peak_summary_counter <- 1L

case_arr <- case_array_integrated

for (unit_key in names(aggregate_indices)) {
  
  parts <- strsplit(unit_key, "\\|\\|")[[1]]
  level_name <- parts[1]
  location_name <- parts[2]
  idx <- aggregate_indices[[unit_key]]
  unit_population <- population_by_unit(idx)
  
  weekly_case_draws <- matrix(
    NA_real_,
    nrow = N_INTEGRATED_DRAWS,
    ncol = 52L
  )
  
  for (wk in 1:52) {
    weekly_case_draws[, wk] <- rowSums(
      case_arr[, idx, wk, drop = FALSE][, , 1, drop = FALSE],
      na.rm = TRUE
    )
  }
  
  peak_week_draw <- max.col(
    weekly_case_draws,
    ties.method = "first"
  )
  
  peak_count_draw <- weekly_case_draws[
    cbind(seq_len(N_INTEGRATED_DRAWS), peak_week_draw)
  ]
  
  peak_incidence_draw <- 100000 * peak_count_draw / unit_population
  
  peak_tab <- as.data.frame(table(peak_week_draw))
  names(peak_tab) <- c("peak_week", "n_draws")
  
  peak_tab$peak_week <- as.integer(as.character(peak_tab$peak_week))
  peak_tab$probability <- peak_tab$n_draws / N_INTEGRATED_DRAWS
  peak_tab$spatial_level <- level_name
  peak_tab$location <- location_name
  
  peak_rows[[peak_counter]] <- tibble::as_tibble(peak_tab) %>%
    dplyr::select(
      .data$spatial_level,
      .data$location,
      .data$peak_week,
      .data$n_draws,
      .data$probability
    )
  
  peak_counter <- peak_counter + 1L
  
  q_peak_inc <- summarise_vector_quantiles(peak_incidence_draw)
  
  peak_summary <- tibble::tibble(
    spatial_level = level_name,
    location = location_name,
    modal_peak_week = as.integer(
      names(which.max(table(peak_week_draw)))
    ),
    probability_modal_peak_week = max(table(peak_week_draw)) /
      N_INTEGRATED_DRAWS
  )
  
  for (nm in FINAL_QNAMES) {
    peak_summary[[paste0("peak_incidence_", nm, "_per100k")]] <-
      q_peak_inc[[nm]]
  }
  
  peak_summary_rows[[peak_summary_counter]] <- peak_summary
  peak_summary_counter <- peak_summary_counter + 1L
}

peak_week_distribution <- dplyr::bind_rows(peak_rows)
peak_summary <- dplyr::bind_rows(peak_summary_rows)

# -----------------------------------------------------------------------------
# 47.13 Long-format quantile submission
# -----------------------------------------------------------------------------

count_quantile_cols <- FINAL_QNAMES

weekly_submission_long <- weekly_submission_wide %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(count_quantile_cols),
    names_to = "quantile_name",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    quantile = dplyr::recode(
      .data$quantile_name,
      q05 = 0.05,
      q10 = 0.10,
      q25 = 0.25,
      median = 0.50,
      q75 = 0.75,
      q90 = 0.90,
      q95 = 0.95
    )
  ) %>%
  dplyr::select(
    .data$model,
    .data$season,
    .data$origin_week,
    .data$target_week,
    .data$horizon,
    .data$outcome,
    .data$spatial_level,
    .data$location,
    .data$population,
    .data$quantile,
    .data$value
  )

EXPECTED_LONG_ROWS <- EXPECTED_WEEKLY_ROWS * 7L

if (nrow(weekly_submission_long) != EXPECTED_LONG_ROWS) {
  stop("Unexpected number of long-format quantile rows.")
}

# -----------------------------------------------------------------------------
# 47.14 Submission quality-control checks
# -----------------------------------------------------------------------------

quantile_matrix <- as.matrix(
  weekly_submission_wide[, FINAL_QNAMES]
)

quantile_order_ok <- all(
  apply(
    quantile_matrix,
    1,
    function(x) all(diff(x) >= 0)
  )
)

submission_qc <- tibble::tibble(
  check = c(
    "weekly_wide_rows",
    "weekly_long_rows",
    "missing_count_quantiles",
    "negative_count_quantiles",
    "quantile_ordering",
    "complete_forecast_weeks",
    "complete_outcomes",
    "number_states",
    "number_zones",
    "number_macro_regions",
    "all_models_converged"
  ),
  value = c(
    as.character(nrow(weekly_submission_wide)),
    as.character(nrow(weekly_submission_long)),
    as.character(sum(is.na(quantile_matrix))),
    as.character(sum(quantile_matrix < 0, na.rm = TRUE)),
    as.character(quantile_order_ok),
    as.character(
      identical(
        sort(unique(weekly_submission_wide$target_week)),
        as.integer(41:52)
      )
    ),
    as.character(
      setequal(
        unique(weekly_submission_wide$outcome),
        c("Cases", "Hospitalizations", "Deaths")
      )
    ),
    as.character(
      dplyr::n_distinct(
        weekly_submission_wide$location[
          weekly_submission_wide$spatial_level == "State"
        ]
      )
    ),
    as.character(
      dplyr::n_distinct(
        weekly_submission_wide$location[
          weekly_submission_wide$spatial_level == "Zone"
        ]
      )
    ),
    as.character(
      dplyr::n_distinct(
        weekly_submission_wide$location[
          weekly_submission_wide$spatial_level == "Macro"
        ]
      )
    ),
    as.character(all(final_diagnostics$convergence_ok))
  )
)

# -----------------------------------------------------------------------------
# 47.15 Parameter summaries for the FINAL origin-40 production fit
# -----------------------------------------------------------------------------

summarise_parameter_final <- function(draw_matrix, variable, model_name) {
  
  x <- draw_matrix[, variable]
  
  tibble::tibble(
    model = model_name,
    parameter = variable,
    mean = mean(x),
    sd = stats::sd(x),
    q05 = unname(stats::quantile(x, 0.05)),
    median = unname(stats::quantile(x, 0.50)),
    q95 = unname(stats::quantile(x, 0.95)),
    probability_positive = mean(x > 0)
  )
}

case_parameter_summary <- dplyr::bind_rows(
  lapply(
    c(
      "beta_memory",
      "beta_national",
      "beta_borrow",
      "beta_phase",
      "beta_interaction"
    ),
    function(v) {
      summarise_parameter_final(
        posterior_matrix,
        v,
        "Cases"
      )
    }
  )
)

hospital_parameter_summary <- dplyr::bind_rows(
  lapply(
    c(
      "beta_case",
      "beta_hosp_memory",
      "beta_phase"
    ),
    function(v) {
      summarise_parameter_final(
        hospital_draws_final,
        v,
        "Hospitalizations"
      )
    }
  )
)

death_parameter_summary <- dplyr::bind_rows(
  lapply(
    c(
      "beta_case",
      "beta_hosp",
      "beta_death_memory",
      "beta_phase"
    ),
    function(v) {
      summarise_parameter_final(
        death_draws_final,
        v,
        "Deaths"
      )
    }
  )
)

final_parameter_summary <- dplyr::bind_rows(
  case_parameter_summary,
  hospital_parameter_summary,
  death_parameter_summary
)

# -----------------------------------------------------------------------------
# 47.16 Forecast metadata
# -----------------------------------------------------------------------------

forecast_metadata <- tibble::tibble(
  model = "PACE-FLU-SEV",
  season = TARGET_SEASON,
  origin_week = FORECAST_ORIGIN,
  first_forecast_week = min(FORECAST_WEEKS),
  last_forecast_week = max(FORECAST_WEEKS),
  forecast_horizon_weeks = length(FORECAST_WEEKS),
  posterior_trajectories = N_INTEGRATED_DRAWS,
  number_states = N_STATES,
  total_population = sum(population_vector),
  growth_scale = growth_scale,
  estimable_connectivity_edges = connectivity$estimable_edges,
  positive_connectivity_edges = connectivity$positive_edges,
  fallback_targets = connectivity$fallback_targets,
  note = paste(
    "Observed through week 40; recursive forecast weeks 41-52.",
    "Seasonal peak/attack-rate targets combine observed weeks 1-40",
    "with simulated weeks 41-52."
  )
)

# -----------------------------------------------------------------------------
# 47.17 Save clean final outputs
# -----------------------------------------------------------------------------

readr::write_csv(
  weekly_submission_wide,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "PACE_FLU_SEV_WEEKLY_SUBMISSION.csv"
  )
)

readr::write_csv(
  weekly_submission_long,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "PACE_FLU_SEV_QUANTILE_SUBMISSION.csv"
  )
)

readr::write_csv(
  direction_probabilities,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "02_DIRECTION_PROBABILITIES.csv"
  )
)

readr::write_csv(
  seasonal_targets,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "03_SEASONAL_TARGETS.csv"
  )
)

readr::write_csv(
  peak_week_distribution,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "04_PEAK_WEEK_DISTRIBUTIONS.csv"
  )
)

readr::write_csv(
  peak_summary,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "05_PEAK_SUMMARY.csv"
  )
)

readr::write_csv(
  final_parameter_summary,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "06_FINAL_PARAMETER_SUMMARY.csv"
  )
)

readr::write_csv(
  forecast_metadata,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "07_FORECAST_METADATA.csv"
  )
)

readr::write_csv(
  submission_qc,
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "08_SUBMISSION_QC.csv"
  )
)

saveRDS(
  list(
    case_array = case_array_integrated,
    hospital_array = hospital_array,
    death_array = death_array,
    state_lookup = state_lookup_final,
    connectivity = W,
    growth_scale = growth_scale,
    diagnostics = final_diagnostics,
    parameters = final_parameter_summary,
    weekly_wide = weekly_submission_wide,
    weekly_long = weekly_submission_long,
    direction_probabilities = direction_probabilities,
    seasonal_targets = seasonal_targets,
    peak_week_distribution = peak_week_distribution,
    peak_summary = peak_summary,
    metadata = forecast_metadata,
    qc = submission_qc
  ),
  file.path(
    INTEGRATED_OUTPUT_DIR,
    "PACE_FLU_SEV_FINAL_RESULTS.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    INTEGRATED_OUTPUT_DIR,
    "09_SESSION_INFO.txt"
  )
)

# -----------------------------------------------------------------------------
# 47.18 Final console report
# -----------------------------------------------------------------------------

cat("\n\n============================================================\n")
cat("PACE-FLU-SEV FINAL SUBMISSION PIPELINE COMPLETE\n")
cat("============================================================\n\n")

cat("MODEL DIAGNOSTICS\n")
print(final_diagnostics, n = Inf)

cat("\nFORECAST METADATA\n")
print(forecast_metadata, n = Inf)

cat("\nSUBMISSION QC\n")
print(submission_qc, n = Inf)

cat("\nFINAL PARAMETER SUMMARY\n")
print(final_parameter_summary, n = Inf)

cat(
  "\nIMPORTANT INTERPRETATION:\n",
  "Peak-week and attack-rate summaries use observed weeks 1-40 plus\n",
  "simulated weeks 41-52. A peak before week 40 is therefore already\n",
  "observed at the forecast origin and is not a prospective peak prediction.\n",
  sep = ""
)

cat("\nOUTPUT DIRECTORY:\n")
cat(
  normalizePath(
    INTEGRATED_OUTPUT_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)




# =============================================================================
# 3. FIGURE 3
#    NORTH VERSUS SOUTH WEEKLY CASE FORECASTS
# =============================================================================

north_south_file <- file.path(
  output_dir,
  "04_NORTH_SOUTH_WEEKLY_FORECASTS.csv"
)

north_south <- read_csv(
  north_south_file,
  show_col_types = FALSE
)


# -----------------------------------------------------------------------------
# Inspect the data
# -----------------------------------------------------------------------------

print(names(north_south))
print(head(north_south))


# -----------------------------------------------------------------------------
# Keep North and South and arrange by week
# -----------------------------------------------------------------------------

north_south_plot_data <- north_south %>%
  filter(
    location %in% c(
      "North",
      "South"
    )
  ) %>%
  arrange(
    location,
    week
  )


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

fig_3 <- ggplot(
  north_south_plot_data,
  aes(
    x = week,
    y = q50,
    group = location
  )
) +
  
  # 90% posterior predictive interval
  geom_ribbon(
    aes(
      ymin = q05,
      ymax = q95,
      fill = location
    ),
    alpha = 0.18,
    colour = NA
  ) +
  
  # Posterior median forecast
  geom_line(
    aes(
      colour = location
    ),
    linewidth = 1.1
  ) +
  
  geom_point(
    aes(
      colour = location
    ),
    size = 2
  ) +
  
  # Observed values for retrospective comparison
  geom_line(
    aes(
      y = actual,
      colour = location
    ),
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  # Forecast origin
  geom_vline(
    xintercept = 40,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  scale_x_continuous(
    breaks = seq(
      min(north_south_plot_data$week),
      max(north_south_plot_data$week),
      by = 1
    )
  ) +
  
  labs(
    title = "North-South Case Forecasts from the Week-40 Origin",
    subtitle = paste0(
      "Solid lines: posterior median; ",
      "shaded regions: 90% predictive intervals; ",
      "dashed lines: observed values"
    ),
    x = "Epidemiological week",
    y = "Weekly cases",
    colour = "Region",
    fill = "Region"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


# Display figure
print(fig_3)


# -----------------------------------------------------------------------------
# Save Figure 3
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_3_North_South_Forecast.png"
  ),
  plot = fig_3,
  width = 9,
  height = 5.5,
  dpi = 300
)


# =============================================================================
# 4. FIGURE 4
#    STATE-LEVEL MEDIAN CASE FORECAST HEATMAP
# =============================================================================

state_file <- file.path(
  output_dir,
  "02_STATE_WEEKLY_FORECASTS.csv"
)

state_forecasts <- read_csv(
  state_file,
  show_col_types = FALSE
)


# -----------------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------------

print(names(state_forecasts))
print(head(state_forecasts))


# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------

state_heatmap_data <- state_forecasts %>%
  filter(
    week >= 41,
    week <= 52
  ) %>%
  mutate(
    state = as.character(state),
    zone = as.character(zone),
    macro_region = as.character(macro_region)
  )


# -----------------------------------------------------------------------------
# Order states geographically rather than alphabetically
#
# States are arranged first by North/South, then geopolitical zone,
# then state name.
# -----------------------------------------------------------------------------

state_order <- state_heatmap_data %>%
  distinct(
    state,
    macro_region,
    zone
  ) %>%
  arrange(
    macro_region,
    zone,
    state
  ) %>%
  pull(state)


state_heatmap_data <- state_heatmap_data %>%
  mutate(
    state = factor(
      state,
      levels = rev(state_order)
    )
  )


# -----------------------------------------------------------------------------
# Plot heatmap
# -----------------------------------------------------------------------------

fig_4 <- ggplot(
  state_heatmap_data,
  aes(
    x = factor(week),
    y = state,
    fill = q50
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.15
  ) +
  
  labs(
    title = "State-Level Median Case Forecasts, Weeks 41-52",
    subtitle = "Week-40 production forecast",
    x = "Forecast week",
    y = "State / FCT",
    fill = "Median\nweekly cases"
  ) +
  
  scale_fill_viridis_c(
    option = "C"
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 7
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  )


# Display figure
print(fig_4)


# -----------------------------------------------------------------------------
# Save Figure 4
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_4_State_Forecast_Heatmap.png"
  ),
  plot = fig_4,
  width = 9,
  height = 10,
  dpi = 300
)


# =============================================================================
# 5. CONFIRM OUTPUTS
# =============================================================================

cat(
  "\n============================================================\n",
  "PACE-FLU-SEV VISUALIZATIONS COMPLETED\n",
  "============================================================\n",
  "\nFigure 3:\n",
  file.path(
    figure_dir,
    "Figure_3_North_South_Forecast.png"
  ),
  "\n\nFigure 4:\n",
  file.path(
    figure_dir,
    "Figure_4_State_Forecast_Heatmap.png"
  ),
  "\n============================================================\n"
)





cat("\n============================================================\n")
cat("END OF FINAL PACE-FLU-SEV SUBMISSION PIPELINE\n")
cat("============================================================\n")


