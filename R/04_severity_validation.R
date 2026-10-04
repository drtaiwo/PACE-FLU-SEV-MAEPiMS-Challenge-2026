# =============================================================================
# PACE-FLU-SEV
# =============================================================================
#
# Phase-Adaptive Connectivity and Epidemic-Memory Forecasting
# with Bayesian Severity Propagation
#
# FINAL INTEGRATED MODEL:
#
#   CASES
#       PACE-FLU M4
#
#            |
#            v
#
#   HOSPITALIZATIONS
#       hierarchical NB severity model
#
#            |
#            v
#
#   DEATHS
#       hierarchical NB mortality model using
#       case burden + hospitalization burden + recent mortality history
#
#
# IMPORTANT:
# ----------
# The cases component is FROZEN M4.
#
# The severity layer does NOT modify the case model.
#
# All transformations are estimated using training data available
# at each forecast origin.
#
# The code performs rolling ONE-WEEK-AHEAD validation because this is
# the operating mode supported by the PACE-FLU validation experiments.
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

SEED <- 20260927

set.seed(SEED)


DATA_FILE <-
  "data/nigeria_flu_weekly_by_state.csv"


OUTPUT_DIR <-
  "PACE_FLU_FINAL_INTEGRATED_SEVERITY"


FIT_DIR <-
  file.path(
    OUTPUT_DIR,
    "fits"
  )


STAN_DIR <-
  file.path(
    OUTPUT_DIR,
    "stan"
  )


FIGURE_DIR <-
  file.path(
    OUTPUT_DIR,
    "figures"
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
  STAN_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


SEASON_LEVELS <-
  c(
    "2023/2024",
    "2024/2025",
    "2025/2026"
  )


# Same definitive origins used in PACE-FLU validation
ORIGINS <-
  c(
    20L,
    24L,
    28L,
    32L,
    36L,
    40L,
    44L,
    48L
  )


# M4 memory
L_MEMORY <- 8L
MEMORY_DECAY <- 0.35


MEMORY_WEIGHTS <-
  exp(
    -MEMORY_DECAY *
      0:(L_MEMORY - 1L)
  )


MEMORY_WEIGHTS <-
  MEMORY_WEIGHTS /
  sum(
    MEMORY_WEIGHTS
  )


MIN_EDGE_N <- 20L


# Stan settings
N_CHAINS <- 4L
N_PARALLEL_CHAINS <- 4L

N_WARMUP <- 1000L
N_SAMPLING <- 1000L

ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15L

INITIAL_STEP_SIZE <- 0.01


# Posterior predictive draws retained
N_PREDICTIVE_DRAWS <- 1000L


# =============================================================================
# 2. GENERAL HELPERS
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


# =============================================================================
# 3. FORECAST SCORING
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
    as.numeric(
      y < lower
    ) +
    
    (2 / alpha) *
    (y - upper) *
    as.numeric(
      y > upper
    )
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


empirical_crps <- function(
    draws,
    y
) {
  
  draws <-
    sort(
      draws[
        is.finite(draws)
      ]
    )
  
  
  n <-
    length(
      draws
    )
  
  
  if (
    n < 2L
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  first_term <-
    mean(
      abs(
        draws - y
      )
    )
  
  
  i <-
    seq_len(n)
  
  
  second_term <-
    sum(
      (
        2 * i -
          n -
          1
      ) *
        draws
    ) /
    n^2
  
  
  first_term -
    second_term
}


summarise_predictive_draws <- function(
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
    
    median = q[4],
    
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
# 4. READ DATA
# =============================================================================

if (
  !file.exists(
    DATA_FILE
  )
) {
  
  stop(
    paste0(
      "Cannot find ",
      DATA_FILE
    )
  )
}


raw_data <-
  readr::read_csv(
    
    DATA_FILE,
    
    show_col_types = FALSE
  )


required_variables <-
  c(
    "season",
    "state",
    "zone",
    "population",
    "epi_week_of_season",
    "cases",
    "hospitalizations",
    "deaths"
  )


missing_variables <-
  setdiff(
    required_variables,
    names(raw_data)
  )


if (
  length(
    missing_variables
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing variables: ",
      paste(
        missing_variables,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 5. PREPARE DATA
# =============================================================================

dat <-
  raw_data %>%
  
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
      as.integer(
        .data$cases
      ),
    
    hospitalizations =
      as.integer(
        .data$hospitalizations
      ),
    
    deaths =
      as.integer(
        .data$deaths
      ),
    
    zone =
      toupper(
        trimws(
          as.character(
            .data$zone
          )
        )
      )
  ) %>%
  
  dplyr::filter(
    
    !is.na(.data$season),
    
    !is.na(.data$state),
    
    !is.na(.data$week),
    
    !is.na(.data$population),
    
    .data$population > 0
  ) %>%
  
  dplyr::mutate(
    
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


# =============================================================================
# 6. STATE INDEX
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


# =============================================================================
# 7. CASE INCIDENCE
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
# 8. NATIONAL CASE SIGNAL
# =============================================================================

national_signal <-
  dat %>%
  
  dplyr::group_by(
    
    .data$season,
    
    .data$season_index,
    
    .data$week
  ) %>%
  
  dplyr::summarise(
    
    national_cases =
      sum(
        .data$cases
      ),
    
    national_population =
      sum(
        .data$population
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
    
    national_signal,
    
    by =
      c(
        "season",
        "season_index",
        "week"
      )
  )


# =============================================================================
# 9. CASE MEMORY
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
      "memory_lag_",
      lag_index
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


memory_variables <-
  paste0(
    "memory_lag_",
    0:(L_MEMORY - 1L)
  )


memory_matrix <-
  as.matrix(
    dat[
      ,
      memory_variables
    ]
  )


dat$memory <-
  as.numeric(
    
    memory_matrix %*%
      MEMORY_WEIGHTS
  )


# =============================================================================
# 10. TEMPORAL FEATURES
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
    
    growth_raw =
      .data$log_incidence -
      dplyr::lag(
        .data$log_incidence,
        1
      ),
    
    cases_next =
      dplyr::lead(
        .data$cases,
        1
      ),
    
    hospitalization_next =
      dplyr::lead(
        .data$hospitalizations,
        1
      ),
    
    deaths_next =
      dplyr::lead(
        .data$deaths,
        1
      ),
    
    log_incidence_next =
      dplyr::lead(
        .data$log_incidence,
        1
      ),
    
    hospitalization_current =
      .data$hospitalizations,
    
    deaths_current =
      .data$deaths,
    
    log_hosp_current =
      log1p(
        .data$hospitalizations
      ),
    
    log_death_current =
      log1p(
        .data$deaths
      ),
    
    recent_hospitalization_rate =
      (
        .data$hospitalizations +
          0.5
      ) /
      (
        .data$cases +
          1
      ),
    
    recent_death_case_rate =
      (
        .data$deaths +
          0.5
      ) /
      (
        .data$cases +
          1
      ),
    
    recent_death_hosp_signal =
      log1p(
        .data$deaths
      ) -
      log1p(
        .data$hospitalizations
      )
  ) %>%
  
  dplyr::ungroup()


# =============================================================================
# 11. CONNECTIVITY ESTIMATION
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
      
      
      coefficient <-
        NA_real_
      
      
      if (
        nrow(edge_data) >=
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
          
          coef_vector <-
            stats::coef(
              edge_fit
            )
          
          
          if (
            "source_current" %in%
            names(coef_vector)
          ) {
            
            coefficient <-
              unname(
                coef_vector[
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
          
          coefficient =
            coefficient
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
        fallback_targets +
        1L
      
      
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
  
  
  diag(W) <- 0
  
  
  list(
    
    W =
      W,
    
    edges =
      edges,
    
    fallback_targets =
      fallback_targets
  )
}


# =============================================================================
# 12. BORROWING
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
      nrow(blocks)
    )
  ) {
    
    indices <-
      which(
        
        data_rows$season ==
          blocks$season[
            block_idx
          ] &
          
          data_rows$week ==
          blocks$week[
            block_idx
          ]
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
# 13. FROZEN CASE MODEL: M4
# =============================================================================

case_stan_code <- '

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

  vector[N] eta;

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


CASE_STAN_FILE <-
  file.path(
    STAN_DIR,
    "PACE_FLU_M4_CASES.stan"
  )


writeLines(
  case_stan_code,
  CASE_STAN_FILE
)


# =============================================================================
# 14. HOSPITALIZATION SEVERITY MODEL
# =============================================================================
#
# H_{i,t+1} ~ NB2(mu_H, phi_H)
#
# log(mu_H) =
#
#   alpha_H
#   + state effect
#   + beta_case_H log(1 + C_{i,t+1})
#   + beta_hosp_H log(1 + H_{it})
#   + beta_phase_H phase_it
#
# The predicted case burden is propagated into hospitalization forecasting.
#
# =============================================================================

hospital_stan_code <- '

data {

  int<lower=1> N;

  int<lower=1> J;

  array[N]
    int<lower=0>
    y;

  array[N]
    int<lower=1,upper=J>
    state_id;

  vector[N]
    log_case_next;

  vector[N]
    log_hosp_current;

  vector[N]
    phase;
}


parameters {

  real alpha;

  vector[J]
    state_raw;

  real<lower=0>
    sigma_state;

  real beta_case;

  real beta_hosp_memory;

  real beta_phase;

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
      -4,
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


  beta_case ~
    normal(
      1,
      0.5
    );


  beta_hosp_memory ~
    normal(
      0,
      1
    );


  beta_phase ~
    normal(
      0,
      1
    );


  log_phi ~
    normal(
      log(5),
      0.75
    );


  phi =
    0.05 +
    exp(
      log_phi
    );


  eta =

      alpha

      + alpha_state[
          state_id
        ]

      + beta_case *
        log_case_next

      + beta_hosp_memory *
        log_hosp_current

      + beta_phase *
        phase;


  y ~
    neg_binomial_2_log(
      eta,
      phi
    );
}
'


HOSPITAL_STAN_FILE <-
  file.path(
    STAN_DIR,
    "PACE_FLU_HOSPITALIZATION_MODEL.stan"
  )


writeLines(
  hospital_stan_code,
  HOSPITAL_STAN_FILE
)


# =============================================================================
# 15. MORTALITY MODEL
# =============================================================================
#
# D_{i,t+1} ~ NB2(mu_D, phi_D)
#
# log(mu_D) =
#
#   alpha_D
#   + state effect
#   + beta_case_D log(1 + C_{i,t+1})
#   + beta_hosp_D log(1 + H_{i,t+1})
#   + beta_death_memory log(1 + D_it)
#   + beta_phase_D phase_it
#
# This does NOT constrain D <= H because the supplied weekly data contain
# some observations where deaths exceed same-week hospitalizations.
#
# =============================================================================

death_stan_code <- '

data {

  int<lower=1> N;

  int<lower=1> J;

  array[N]
    int<lower=0>
    y;

  array[N]
    int<lower=1,upper=J>
    state_id;

  vector[N]
    log_case_next;

  vector[N]
    log_hosp_next;

  vector[N]
    log_death_current;

  vector[N]
    phase;
}


parameters {

  real alpha;

  vector[J]
    state_raw;

  real<lower=0>
    sigma_state;

  real beta_case;

  real beta_hosp;

  real beta_death_memory;

  real beta_phase;

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
      -5,
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


  beta_case ~
    normal(
      0.5,
      0.5
    );


  beta_hosp ~
    normal(
      0.5,
      0.5
    );


  beta_death_memory ~
    normal(
      0,
      1
    );


  beta_phase ~
    normal(
      0,
      1
    );


  log_phi ~
    normal(
      log(5),
      0.75
    );


  phi =
    0.05 +
    exp(
      log_phi
    );


  eta =

      alpha

      + alpha_state[
          state_id
        ]

      + beta_case *
        log_case_next

      + beta_hosp *
        log_hosp_next

      + beta_death_memory *
        log_death_current

      + beta_phase *
        phase;


  y ~
    neg_binomial_2_log(
      eta,
      phi
    );
}
'


DEATH_STAN_FILE <-
  file.path(
    STAN_DIR,
    "PACE_FLU_DEATH_MODEL.stan"
  )


writeLines(
  death_stan_code,
  DEATH_STAN_FILE
)


# =============================================================================
# 16. COMPILE MODELS
# =============================================================================

cat(
  "\nCompiling case model...\n"
)


case_model <-
  cmdstanr::cmdstan_model(
    CASE_STAN_FILE
  )


cat(
  "\nCompiling hospitalization model...\n"
)


hospital_model <-
  cmdstanr::cmdstan_model(
    HOSPITAL_STAN_FILE
  )


cat(
  "\nCompiling mortality model...\n"
)


death_model <-
  cmdstanr::cmdstan_model(
    DEATH_STAN_FILE
  )


# =============================================================================
# 17. STORAGE
# =============================================================================

all_forecasts <-
  list()


all_diagnostics <-
  list()


all_parameters <-
  list()


forecast_counter <- 1L
diagnostic_counter <- 1L
parameter_counter <- 1L


# =============================================================================
# 18. MAIN VALIDATION LOOP
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
      "ORIGIN: ",
      origin,
      "\n",
      sep = ""
    )
    
    cat(
      "TARGET WEEK: ",
      origin + 1L,
      "\n",
      sep = ""
    )
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    
    # =========================================================================
    # 18.1 TRAINING SET
    # =========================================================================
    
    train <-
      dat %>%
      
      dplyr::filter(
        
        !is.na(
          .data$cases_next
        ),
        
        !is.na(
          .data$hospitalization_next
        ),
        
        !is.na(
          .data$deaths_next
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
    
    
    if (
      nrow(train) <
      200
    ) {
      
      warning(
        paste0(
          "Insufficient training data for ",
          target_season,
          " origin ",
          origin
        )
      )
      
      next
    }
    
    
    # =========================================================================
    # 18.2 TRAIN-ONLY PHASE SCALE
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
    # 18.3 CONNECTIVITY
    # =========================================================================
    
    connectivity <-
      estimate_connectivity(
        
        train,
        
        N_STATES
      )
    
    
    W <-
      connectivity$W
    
    
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
          
          .data$hospitalization_next,
          
          .data$deaths_next,
          
          .data$memory,
          
          .data$national_log_incidence,
          
          .data$borrowing,
          
          .data$p_growth,
          
          .data$borrow_phase,
          
          .data$log_hosp_current,
          
          .data$log_death_current
        )
      )
    
    
    # =========================================================================
    # 18.4 FORECAST-ORIGIN ROWS
    # =========================================================================
    
    origin_data <-
      dat %>%
      
      dplyr::filter(
        
        .data$season ==
          target_season,
        
        .data$week ==
          origin
      ) %>%
      
      dplyr::arrange(
        .data$state_id
      )
    
    
    target_data <-
      dat %>%
      
      dplyr::filter(
        
        .data$season ==
          target_season,
        
        .data$week ==
          origin + 1L
      ) %>%
      
      dplyr::arrange(
        .data$state_id
      )
    
    
    if (
      nrow(origin_data) !=
      N_STATES ||
      nrow(target_data) !=
      N_STATES
    ) {
      
      warning(
        "Incomplete forecast origin/target."
      )
      
      next
    }
    
    
    # =========================================================================
    # 18.5 ORIGIN PHASE
    # =========================================================================
    
    origin_data <-
      origin_data %>%
      
      dplyr::mutate(
        
        p_growth =
          plogis(
            .data$growth_raw /
              growth_scale
          )
      )
    
    
    # =========================================================================
    # 18.6 ORIGIN BORROWING
    # =========================================================================
    
    current_log_incidence <-
      origin_data$log_incidence
    
    
    origin_borrowing <-
      as.numeric(
        
        W %*%
          current_log_incidence
      )
    
    
    origin_data$borrowing <-
      origin_borrowing
    
    
    origin_data$borrow_phase <-
      origin_data$borrowing *
      origin_data$p_growth
    
    
    # =========================================================================
    # 18.7 CASE MODEL DATA
    # =========================================================================
    
    case_data <-
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
    # 18.8 HOSPITAL MODEL DATA
    #
    # Training is allowed to use observed cases_next because they are outcomes
    # from historical training rows.
    #
    # At forecast time, cases_next will be replaced by posterior predicted
    # cases, propagating case uncertainty.
    # =========================================================================
    
    hospital_data <-
      list(
        
        N =
          nrow(train),
        
        J =
          N_STATES,
        
        y =
          as.integer(
            train$hospitalization_next
          ),
        
        state_id =
          as.integer(
            train$state_id
          ),
        
        log_case_next =
          log1p(
            train$cases_next
          ),
        
        log_hosp_current =
          as.vector(
            train$log_hosp_current
          ),
        
        phase =
          as.vector(
            train$p_growth
          )
      )
    
    
    # =========================================================================
    # 18.9 DEATH MODEL DATA
    # =========================================================================
    
    death_data <-
      list(
        
        N =
          nrow(train),
        
        J =
          N_STATES,
        
        y =
          as.integer(
            train$deaths_next
          ),
        
        state_id =
          as.integer(
            train$state_id
          ),
        
        log_case_next =
          log1p(
            train$cases_next
          ),
        
        log_hosp_next =
          log1p(
            train$hospitalization_next
          ),
        
        log_death_current =
          as.vector(
            train$log_death_current
          ),
        
        phase =
          as.vector(
            train$p_growth
          )
      )
    
    
    # =========================================================================
    # 18.10 CHECKPOINT
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
          "integrated_",
          season_tag,
          "_origin_",
          origin,
          ".rds"
        )
      )
    
    
    # =========================================================================
    # 18.11 FIT MODELS
    # =========================================================================
    
    if (
      file.exists(
        checkpoint_file
      )
    ) {
      
      cat(
        "Loading checkpoint...\n"
      )
      
      
      saved_fit <-
        readRDS(
          checkpoint_file
        )
      
      
      case_draws <-
        saved_fit$case_draws
      
      
      hospital_draws <-
        saved_fit$hospital_draws
      
      
      death_draws <-
        saved_fit$death_draws
      
      
      diagnostics_current <-
        saved_fit$diagnostics
      
      
    } else {
      
      # =======================================================================
      # CASE MODEL
      # =======================================================================
      
      cat(
        "Fitting frozen PACE-FLU M4 case model...\n"
      )
      
      
      fit_case <-
        case_model$sample(
          
          data =
            case_data,
          
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
          
          adapt_delta =
            ADAPT_DELTA,
          
          max_treedepth =
            MAX_TREEDEPTH,
          
          step_size =
            INITIAL_STEP_SIZE,
          
          init =
            0,
          
          refresh =
            500
        )
      
      
      # =======================================================================
      # HOSPITAL MODEL
      # =======================================================================
      
      cat(
        "Fitting hospitalization severity model...\n"
      )
      
      
      fit_hospital <-
        hospital_model$sample(
          
          data =
            hospital_data,
          
          seed =
            SEED +
            100000L +
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
          
          adapt_delta =
            ADAPT_DELTA,
          
          max_treedepth =
            MAX_TREEDEPTH,
          
          step_size =
            INITIAL_STEP_SIZE,
          
          init =
            0,
          
          refresh =
            500
        )
      
      
      # =======================================================================
      # DEATH MODEL
      # =======================================================================
      
      cat(
        "Fitting mortality model...\n"
      )
      
      
      fit_death <-
        death_model$sample(
          
          data =
            death_data,
          
          seed =
            SEED +
            200000L +
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
          
          adapt_delta =
            ADAPT_DELTA,
          
          max_treedepth =
            MAX_TREEDEPTH,
          
          step_size =
            INITIAL_STEP_SIZE,
          
          init =
            0,
          
          refresh =
            500
        )
      
      
      # =======================================================================
      # DIAGNOSTICS
      # =======================================================================
      
      get_diagnostics <- function(
    fit,
    model_name
      ) {
        
        sm <-
          fit$summary()
        
        
        ds <-
          fit$diagnostic_summary()
        
        
        tibble::tibble(
          
          season =
            target_season,
          
          origin =
            origin,
          
          model =
            model_name,
          
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
              ds$num_divergent,
              na.rm = TRUE
            ),
          
          treedepth_hits =
            sum(
              ds$num_max_treedepth,
              na.rm = TRUE
            )
        )
      }
      
      
      diagnostics_current <-
        dplyr::bind_rows(
          
          get_diagnostics(
            fit_case,
            "Cases"
          ),
          
          get_diagnostics(
            fit_hospital,
            "Hospitalizations"
          ),
          
          get_diagnostics(
            fit_death,
            "Deaths"
          )
        )
      
      
      # =======================================================================
      # POSTERIOR DRAWS
      # =======================================================================
      
      case_draws <-
        fit_case$draws(
          format = "matrix"
        )
      
      
      hospital_draws <-
        fit_hospital$draws(
          format = "matrix"
        )
      
      
      death_draws <-
        fit_death$draws(
          format = "matrix"
        )
      
      
      # =======================================================================
      # CHECKPOINT
      # =======================================================================
      
      saveRDS(
        
        list(
          
          case_draws =
            case_draws,
          
          hospital_draws =
            hospital_draws,
          
          death_draws =
            death_draws,
          
          diagnostics =
            diagnostics_current
        ),
        
        checkpoint_file
      )
      
      
      rm(
        fit_case,
        fit_hospital,
        fit_death
      )
      
      
      gc()
    }
    
    
    # =========================================================================
    # 18.12 SAVE DIAGNOSTICS
    # =========================================================================
    
    all_diagnostics[[
      diagnostic_counter
    ]] <-
      diagnostics_current
    
    
    diagnostic_counter <-
      diagnostic_counter + 1L
    
    
    # =========================================================================
    # 18.13 MATCH POSTERIOR DRAWS
    # =========================================================================
    
    n_available <-
      min(
        nrow(case_draws),
        nrow(hospital_draws),
        nrow(death_draws)
      )
    
    
    n_use <-
      min(
        N_PREDICTIVE_DRAWS,
        n_available
      )
    
    
    draw_indices <-
      unique(
        round(
          seq(
            1,
            n_available,
            length.out =
              n_use
          )
        )
      )
    
    
    case_draws <-
      case_draws[
        draw_indices,
        ,
        drop = FALSE
      ]
    
    
    hospital_draws <-
      hospital_draws[
        draw_indices,
        ,
        drop = FALSE
      ]
    
    
    death_draws <-
      death_draws[
        draw_indices,
        ,
        drop = FALSE
      ]
    
    
    N_DRAWS <-
      length(
        draw_indices
      )
    
    
    # =========================================================================
    # 18.14 STORAGE FOR PREDICTIONS
    # =========================================================================
    
    case_predictions <-
      matrix(
        
        NA_real_,
        
        nrow =
          N_DRAWS,
        
        ncol =
          N_STATES
      )
    
    
    hospital_predictions <-
      matrix(
        
        NA_real_,
        
        nrow =
          N_DRAWS,
        
        ncol =
          N_STATES
      )
    
    
    death_predictions <-
      matrix(
        
        NA_real_,
        
        nrow =
          N_DRAWS,
        
        ncol =
          N_STATES
      )
    
    
    # =========================================================================
    # 18.15 POSTERIOR PREDICTIVE CASCADE
    # =========================================================================
    
    for (
      draw_idx in
      seq_len(
        N_DRAWS
      )
    ) {
      
      # =======================================================================
      # CASE PARAMETERS
      # =======================================================================
      
      case_alpha <-
        case_draws[
          draw_idx,
          "alpha"
        ]
      
      
      case_state_effect <-
        sapply(
          
          seq_len(
            N_STATES
          ),
          
          function(j) {
            
            case_draws[
              draw_idx,
              paste0(
                "alpha_state[",
                j,
                "]"
              )
            ]
          }
        )
      
      
      case_phi <-
        0.05 +
        exp(
          case_draws[
            draw_idx,
            "log_phi"
          ]
        )
      
      
      # =======================================================================
      # CASE EXPECTATION
      # =======================================================================
      
      eta_case <-
        
        log(
          population_vector /
            100000
        ) +
        
        case_alpha +
        
        case_state_effect +
        
        case_draws[
          draw_idx,
          "beta_memory"
        ] *
        origin_data$memory +
        
        case_draws[
          draw_idx,
          "beta_national"
        ] *
        origin_data$national_log_incidence +
        
        case_draws[
          draw_idx,
          "beta_borrow"
        ] *
        origin_data$borrowing +
        
        case_draws[
          draw_idx,
          "beta_phase"
        ] *
        origin_data$p_growth +
        
        case_draws[
          draw_idx,
          "beta_interaction"
        ] *
        origin_data$borrow_phase
      
      
      mu_case <-
        exp(
          eta_case
        )
      
      
      predicted_cases <-
        stats::rnbinom(
          
          N_STATES,
          
          mu =
            mu_case,
          
          size =
            case_phi
        )
      
      
      case_predictions[
        draw_idx,
      ] <-
        predicted_cases
      
      
      # =======================================================================
      # HOSPITALIZATION PARAMETERS
      # =======================================================================
      
      hosp_alpha <-
        hospital_draws[
          draw_idx,
          "alpha"
        ]
      
      
      hosp_state_effect <-
        sapply(
          
          seq_len(
            N_STATES
          ),
          
          function(j) {
            
            hospital_draws[
              draw_idx,
              paste0(
                "alpha_state[",
                j,
                "]"
              )
            ]
          }
        )
      
      
      hosp_phi <-
        0.05 +
        exp(
          hospital_draws[
            draw_idx,
            "log_phi"
          ]
        )
      
      
      # =======================================================================
      # HOSPITALIZATION EXPECTATION
      #
      # Notice log1p(predicted_cases):
      # uncertainty from the case model is propagated.
      # =======================================================================
      
      eta_hosp <-
        
        hosp_alpha +
        
        hosp_state_effect +
        
        hospital_draws[
          draw_idx,
          "beta_case"
        ] *
        log1p(
          predicted_cases
        ) +
        
        hospital_draws[
          draw_idx,
          "beta_hosp_memory"
        ] *
        origin_data$log_hosp_current +
        
        hospital_draws[
          draw_idx,
          "beta_phase"
        ] *
        origin_data$p_growth
      
      
      mu_hosp <-
        exp(
          eta_hosp
        )
      
      
      predicted_hospitalizations <-
        stats::rnbinom(
          
          N_STATES,
          
          mu =
            mu_hosp,
          
          size =
            hosp_phi
        )
      
      
      hospital_predictions[
        draw_idx,
      ] <-
        predicted_hospitalizations
      
      
      # =======================================================================
      # DEATH PARAMETERS
      # =======================================================================
      
      death_alpha <-
        death_draws[
          draw_idx,
          "alpha"
        ]
      
      
      death_state_effect <-
        sapply(
          
          seq_len(
            N_STATES
          ),
          
          function(j) {
            
            death_draws[
              draw_idx,
              paste0(
                "alpha_state[",
                j,
                "]"
              )
            ]
          }
        )
      
      
      death_phi <-
        0.05 +
        exp(
          death_draws[
            draw_idx,
            "log_phi"
          ]
        )
      
      
      # =======================================================================
      # DEATH EXPECTATION
      #
      # Both case and hospitalization uncertainty propagate into mortality.
      # =======================================================================
      
      eta_death <-
        
        death_alpha +
        
        death_state_effect +
        
        death_draws[
          draw_idx,
          "beta_case"
        ] *
        log1p(
          predicted_cases
        ) +
        
        death_draws[
          draw_idx,
          "beta_hosp"
        ] *
        log1p(
          predicted_hospitalizations
        ) +
        
        death_draws[
          draw_idx,
          "beta_death_memory"
        ] *
        origin_data$log_death_current +
        
        death_draws[
          draw_idx,
          "beta_phase"
        ] *
        origin_data$p_growth
      
      
      mu_death <-
        exp(
          eta_death
        )
      
      
      predicted_deaths <-
        stats::rnbinom(
          
          N_STATES,
          
          mu =
            mu_death,
          
          size =
            death_phi
        )
      
      
      death_predictions[
        draw_idx,
      ] <-
        predicted_deaths
    }
    
    
    # =========================================================================
    # 18.16 STATE-LEVEL SCORING
    # =========================================================================
    
    for (
      state_idx in
      seq_len(
        N_STATES
      )
    ) {
      
      common_info <-
        tibble::tibble(
          
          season =
            target_season,
          
          origin =
            origin,
          
          target_week =
            origin + 1L,
          
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
        )
      
      
      # Cases
      all_forecasts[[
        forecast_counter
      ]] <-
        
        common_info %>%
        
        dplyr::mutate(
          outcome = "Cases"
        ) %>%
        
        dplyr::bind_cols(
          
          summarise_predictive_draws(
            
            case_predictions[
              ,
              state_idx
            ],
            
            target_data$cases[
              state_idx
            ]
          )
        )
      
      
      forecast_counter <-
        forecast_counter + 1L
      
      
      # Hospitalizations
      all_forecasts[[
        forecast_counter
      ]] <-
        
        common_info %>%
        
        dplyr::mutate(
          outcome =
            "Hospitalizations"
        ) %>%
        
        dplyr::bind_cols(
          
          summarise_predictive_draws(
            
            hospital_predictions[
              ,
              state_idx
            ],
            
            target_data$hospitalizations[
              state_idx
            ]
          )
        )
      
      
      forecast_counter <-
        forecast_counter + 1L
      
      
      # Deaths
      all_forecasts[[
        forecast_counter
      ]] <-
        
        common_info %>%
        
        dplyr::mutate(
          outcome =
            "Deaths"
        ) %>%
        
        dplyr::bind_cols(
          
          summarise_predictive_draws(
            
            death_predictions[
              ,
              state_idx
            ],
            
            target_data$deaths[
              state_idx
            ]
          )
        )
      
      
      forecast_counter <-
        forecast_counter + 1L
    }
    
    
    # =========================================================================
    # 18.17 PARAMETER SUMMARIES
    # =========================================================================
    
    summarize_parameter <- function(
    draws,
    variable,
    model_name
    ) {
      
      x <-
        draws[
          ,
          variable
        ]
      
      
      tibble::tibble(
        
        season =
          target_season,
        
        origin =
          origin,
        
        model =
          model_name,
        
        parameter =
          variable,
        
        mean =
          mean(x),
        
        sd =
          stats::sd(x),
        
        q05 =
          unname(
            stats::quantile(
              x,
              0.05
            )
          ),
        
        q50 =
          unname(
            stats::quantile(
              x,
              0.50
            )
          ),
        
        q95 =
          unname(
            stats::quantile(
              x,
              0.95
            )
          ),
        
        probability_positive =
          mean(
            x > 0
          )
      )
    }
    
    
    current_parameters <-
      dplyr::bind_rows(
        
        summarize_parameter(
          case_draws,
          "beta_memory",
          "Cases"
        ),
        
        summarize_parameter(
          case_draws,
          "beta_national",
          "Cases"
        ),
        
        summarize_parameter(
          case_draws,
          "beta_borrow",
          "Cases"
        ),
        
        summarize_parameter(
          case_draws,
          "beta_phase",
          "Cases"
        ),
        
        summarize_parameter(
          case_draws,
          "beta_interaction",
          "Cases"
        ),
        
        summarize_parameter(
          hospital_draws,
          "beta_case",
          "Hospitalizations"
        ),
        
        summarize_parameter(
          hospital_draws,
          "beta_hosp_memory",
          "Hospitalizations"
        ),
        
        summarize_parameter(
          hospital_draws,
          "beta_phase",
          "Hospitalizations"
        ),
        
        summarize_parameter(
          death_draws,
          "beta_case",
          "Deaths"
        ),
        
        summarize_parameter(
          death_draws,
          "beta_hosp",
          "Deaths"
        ),
        
        summarize_parameter(
          death_draws,
          "beta_death_memory",
          "Deaths"
        ),
        
        summarize_parameter(
          death_draws,
          "beta_phase",
          "Deaths"
        )
      )
    
    
    all_parameters[[
      parameter_counter
    ]] <-
      current_parameters
    
    
    parameter_counter <-
      parameter_counter + 1L
    
    
    # =========================================================================
    # MEMORY CLEANUP
    # =========================================================================
    
    rm(
      case_draws,
      hospital_draws,
      death_draws,
      case_predictions,
      hospital_predictions,
      death_predictions
    )
    
    
    gc()
  }
}


# =============================================================================
# 19. COMBINE RESULTS
# =============================================================================

forecast_results <-
  dplyr::bind_rows(
    all_forecasts
  )


diagnostic_results <-
  dplyr::bind_rows(
    all_diagnostics
  )


parameter_results <-
  dplyr::bind_rows(
    all_parameters
  )


# =============================================================================
# 20. OVERALL PERFORMANCE
# =============================================================================

overall_performance <-
  forecast_results %>%
  
  dplyr::group_by(
    .data$outcome
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
    
    .groups =
      "drop"
  )


# =============================================================================
# 21. PERFORMANCE BY SEASON
# =============================================================================

performance_by_season <-
  forecast_results %>%
  
  dplyr::group_by(
    
    .data$outcome,
    
    .data$season
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


# =============================================================================
# 22. PERFORMANCE BY ORIGIN
# =============================================================================

performance_by_origin <-
  forecast_results %>%
  
  dplyr::group_by(
    
    .data$outcome,
    
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


# =============================================================================
# 23. AGGREGATED NATIONAL FORECASTS
# =============================================================================
#
# IMPORTANT:
#
# State quantiles cannot simply be summed.
#
# For definitive aggregate uncertainty we would retain and aggregate the
# posterior draws themselves.
#
# This validation file therefore reports state-level scoring as the primary
# evaluation. National summaries below are based on observed/median sums and
# should be treated as descriptive point-forecast summaries.
#
# =============================================================================

national_point_forecasts <-
  forecast_results %>%
  
  dplyr::group_by(
    
    .data$outcome,
    
    .data$season,
    
    .data$origin,
    
    .data$target_week
  ) %>%
  
  dplyr::summarise(
    
    actual =
      sum(
        .data$actual
      ),
    
    predicted_median_sum =
      sum(
        .data$median
      ),
    
    error =
      .data$predicted_median_sum -
      .data$actual,
    
    absolute_error =
      abs(
        .data$error
      ),
    
    .groups =
      "drop"
  )


# =============================================================================
# 24. SAVE OUTPUTS
# =============================================================================

readr::write_csv(
  
  forecast_results,
  
  file.path(
    OUTPUT_DIR,
    "01_STATE_LEVEL_FORECASTS.csv"
  )
)


readr::write_csv(
  
  overall_performance,
  
  file.path(
    OUTPUT_DIR,
    "02_OVERALL_PERFORMANCE.csv"
  )
)


readr::write_csv(
  
  performance_by_season,
  
  file.path(
    OUTPUT_DIR,
    "03_PERFORMANCE_BY_SEASON.csv"
  )
)


readr::write_csv(
  
  performance_by_origin,
  
  file.path(
    OUTPUT_DIR,
    "04_PERFORMANCE_BY_ORIGIN.csv"
  )
)


readr::write_csv(
  
  diagnostic_results,
  
  file.path(
    OUTPUT_DIR,
    "05_MODEL_DIAGNOSTICS.csv"
  )
)


readr::write_csv(
  
  parameter_results,
  
  file.path(
    OUTPUT_DIR,
    "06_PARAMETER_SUMMARIES.csv"
  )
)


readr::write_csv(
  
  national_point_forecasts,
  
  file.path(
    OUTPUT_DIR,
    "07_NATIONAL_POINT_FORECASTS.csv"
  )
)


# =============================================================================
# 25. SEVERITY PARAMETER SUMMARY ACROSS ORIGINS
# =============================================================================

parameter_overall <-
  parameter_results %>%
  
  dplyr::group_by(
    
    .data$model,
    
    .data$parameter
  ) %>%
  
  dplyr::summarise(
    
    origins =
      dplyr::n(),
    
    mean_posterior_median =
      mean(
        .data$q50
      ),
    
    median_posterior_median =
      stats::median(
        .data$q50
      ),
    
    minimum_posterior_median =
      min(
        .data$q50
      ),
    
    maximum_posterior_median =
      max(
        .data$q50
      ),
    
    proportion_positive =
      mean(
        .data$q50 > 0
      ),
    
    proportion_90CrI_above_zero =
      mean(
        .data$q05 > 0
      ),
    
    mean_probability_positive =
      mean(
        .data$probability_positive
      ),
    
    .groups =
      "drop"
  )


readr::write_csv(
  
  parameter_overall,
  
  file.path(
    OUTPUT_DIR,
    "08_PARAMETER_EVIDENCE_ACROSS_ORIGINS.csv"
  )
)


# =============================================================================
# 26. FIGURE: PERFORMANCE BY OUTCOME
# =============================================================================

p_mae <-
  ggplot2::ggplot(
    
    performance_by_origin,
    
    ggplot2::aes(
      x = origin,
      y = MAE,
      group = season,
      linetype = season
    )
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::facet_wrap(
    
    ~ outcome,
    
    scales =
      "free_y"
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU-SEV Rolling Forecast Accuracy",
    
    subtitle =
      "One-week-ahead validation",
    
    x =
      "Forecast origin",
    
    y =
      "Mean absolute error",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "01_MAE_BY_ORIGIN.png"
  ),
  
  p_mae,
  
  width = 10,
  
  height = 7,
  
  dpi = 300
)


# =============================================================================
# 27. FIGURE: COVERAGE
# =============================================================================

p_coverage <-
  performance_by_origin %>%
  
  ggplot2::ggplot(
    
    ggplot2::aes(
      x = origin,
      y = coverage90,
      group = season,
      linetype = season
    )
  ) +
  
  ggplot2::geom_hline(
    
    yintercept =
      0.90,
    
    linetype =
      2
  ) +
  
  ggplot2::geom_line() +
  
  ggplot2::geom_point() +
  
  ggplot2::facet_wrap(
    ~ outcome
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
      "Empirical 90% Predictive-Interval Coverage",
    
    x =
      "Forecast origin",
    
    y =
      "Coverage",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "02_COVERAGE90_BY_ORIGIN.png"
  ),
  
  p_coverage,
  
  width = 10,
  
  height = 7,
  
  dpi = 300
)


# =============================================================================
# 28. FIGURE: BIAS
# =============================================================================

p_bias <-
  performance_by_origin %>%
  
  ggplot2::ggplot(
    
    ggplot2::aes(
      x = origin,
      y = bias,
      group = season,
      linetype = season
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
  
  ggplot2::facet_wrap(
    
    ~ outcome,
    
    scales =
      "free_y"
  ) +
  
  ggplot2::labs(
    
    title =
      "PACE-FLU-SEV Forecast Bias",
    
    subtitle =
      "Positive values indicate overprediction",
    
    x =
      "Forecast origin",
    
    y =
      "Mean forecast minus observation",
    
    linetype =
      "Season"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  
  file.path(
    FIGURE_DIR,
    "03_BIAS_BY_ORIGIN.png"
  ),
  
  p_bias,
  
  width = 10,
  
  height = 7,
  
  dpi = 300
)


# =============================================================================
# 29. SAVE COMPLETE WORKSPACE OBJECT
# =============================================================================

saveRDS(
  
  list(
    
    state_forecasts =
      forecast_results,
    
    overall_performance =
      overall_performance,
    
    performance_by_season =
      performance_by_season,
    
    performance_by_origin =
      performance_by_origin,
    
    diagnostics =
      diagnostic_results,
    
    parameter_summaries =
      parameter_results,
    
    parameter_evidence =
      parameter_overall,
    
    national_point_forecasts =
      national_point_forecasts
  ),
  
  file.path(
    OUTPUT_DIR,
    "PACE_FLU_SEV_COMPLETE_RESULTS.rds"
  )
)


# =============================================================================
# 30. SESSION INFORMATION
# =============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      OUTPUT_DIR,
      "09_SESSION_INFO.txt"
    )
)


# =============================================================================
# 31. PRINT FINAL RESULTS
# =============================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "PACE-FLU-SEV VALIDATION COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nOVERALL PERFORMANCE\n\n"
)


print(
  overall_performance,
  n = Inf
)


cat(
  "\n\nPERFORMANCE BY SEASON\n\n"
)


print(
  performance_by_season,
  n = Inf
)


cat(
  "\n\nPARAMETER EVIDENCE\n\n"
)


print(
  parameter_overall,
  n = Inf
)


cat(
  "\n\nMODEL DIAGNOSTICS\n\n"
)


print(
  diagnostic_results,
  n = Inf
)


cat(
  "\n\nOUTPUT DIRECTORY:\n"
)


cat(
  normalizePath(
    
    OUTPUT_DIR,
    
    winslash = "/",
    
    mustWork = FALSE
  )
)


cat(
  "\n\n============================================================\n"
)

cat(
  "END OF PACE-FLU-SEV\n"
)

cat(
  "============================================================\n"
)