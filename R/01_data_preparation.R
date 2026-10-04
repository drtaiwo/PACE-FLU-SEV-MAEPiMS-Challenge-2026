# =============================================================================
# PACE-FLU: DATA PREPARATION AND QUALITY CONTROL
# =============================================================================
# Uses only the supplied MAEPiMS synthetic state-level surveillance dataset.
# Expected path: data/nigeria_flu_weekly_by_state.csv
# This script performs reproducible checks and writes a cleaned analysis file.
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)

DATA_FILE <- "data/nigeria_flu_weekly_by_state.csv"
OUTPUT_DIR <- "PACE_FLU_DATA_PREPARATION"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(DATA_FILE)) stop("Missing data file: ", DATA_FILE)

normalize_season <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("2023/24", "2023/2024") ~ "2023/2024",
    x %in% c("2024/25", "2024/2025") ~ "2024/2025",
    x %in% c("2025/26", "2025/2026") ~ "2025/2026",
    TRUE ~ NA_character_
  )
}

raw <- readr::read_csv(DATA_FILE, show_col_types = FALSE)
required <- c("season","state","zone","population","epi_week_of_season",
              "cases","hospitalizations","deaths")
missing <- setdiff(required, names(raw))
if (length(missing) > 0L) stop("Missing required columns: ", paste(missing, collapse=", "))

dat <- raw %>%
  mutate(
    season = normalize_season(season),
    state = as.character(state),
    zone = toupper(trimws(as.character(zone))),
    population = as.numeric(population),
    epi_week_of_season = as.integer(epi_week_of_season),
    cases = as.numeric(cases),
    hospitalizations = as.numeric(hospitalizations),
    deaths = as.numeric(deaths),
    cases_per_100k_calc = 100000 * cases / population,
    hosp_per_100k_calc = 100000 * hospitalizations / population,
    deaths_per_100k_calc = 100000 * deaths / population,
    macro_region = case_when(
      zone %in% c("NW","NE","NC") ~ "North",
      zone %in% c("SW","SE","SS") ~ "South",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(season, state, epi_week_of_season)

if (anyNA(dat$season)) stop("Unrecognized season label detected.")
if (any(!is.finite(dat$population) | dat$population <= 0)) stop("Invalid population values.")
if (any(dat$cases < 0 | dat$hospitalizations < 0 | dat$deaths < 0, na.rm=TRUE)) stop("Negative outcomes detected.")

duplicates <- dat %>% count(season, state, epi_week_of_season) %>% filter(n > 1)
if (nrow(duplicates) > 0L) stop("Duplicate state-season-week rows detected.")

structure_summary <- dat %>%
  count(season, state, name="n_weeks") %>%
  summarise(n_state_seasons=n(), min_weeks=min(n_weeks), max_weeks=max(n_weeks))

qc <- tibble::tibble(
  rows = nrow(dat),
  states = n_distinct(dat$state),
  seasons = n_distinct(dat$season),
  min_week = min(dat$epi_week_of_season),
  max_week = max(dat$epi_week_of_season),
  missing_cases = sum(is.na(dat$cases)),
  missing_hospitalizations = sum(is.na(dat$hospitalizations)),
  missing_deaths = sum(is.na(dat$deaths)),
  hospitalizations_gt_cases = sum(dat$hospitalizations > dat$cases, na.rm=TRUE),
  deaths_gt_hospitalizations = sum(dat$deaths > dat$hospitalizations, na.rm=TRUE)
)

readr::write_csv(dat, file.path(OUTPUT_DIR, "01_clean_state_data.csv"))
readr::write_csv(qc, file.path(OUTPUT_DIR, "02_quality_control_summary.csv"))
readr::write_csv(structure_summary, file.path(OUTPUT_DIR, "03_structure_summary.csv"))
capture.output(sessionInfo(), file=file.path(OUTPUT_DIR, "04_session_info.txt"))
print(qc)
print(structure_summary)
cat("\nData preparation complete.\n")
