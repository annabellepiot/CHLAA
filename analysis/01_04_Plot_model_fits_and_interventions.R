## 01_04  Publication figure — model fits + interventions
## Produces a composite figure with a large Kirotshe panel (with full
## intervention labels and a shared legend) and smaller facets for every
## other health zone.
##
## Run interactively or via:
##   conda activate Renv
##   Rscript analysis/01_04_Plot_model_fits_and_interventions.R

library(chlaa)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# ── paths ────────────────────────────────────────────────────────────────
base_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA"
data_dir <- file.path(base_dir, "analysis", "data")
rds_dir <- file.path(base_dir, "figures", ".rds files")
fig_dir <- file.path(base_dir, "figures")

# ── colour palette for intervention types ────────────────────────────────
# Bright, high-contrast colours – each maximally distinguishable
intervention_colours <- c(
  "CTC"           = "#00BFC4",
  "ORC"           = "#f86d6d",
  "CATI"          = "#7B68EE",
  "Hygiene"       = "#FF69B4",
  "Chlorination"  = "#32CD32",
  "Latrines"      = "#FFD700",
  "Vaccination"   = "#FF8C00"
)

# nice labels for the legend
intervention_labels <- names(intervention_colours)

# ── load intervention parameters ─────────────────────────────────────────
hz_params <- read_csv(file.path(data_dir, "hz_parameters.csv"),
  show_col_types = FALSE
)

# pivot to one row per HZ with start/end dates for each intervention
parse_interventions <- function(hz_params) {
  # intervention type mapping  (prefix in parameter name → display label)
  type_map <- c(
    ctc   = "CTC",
    orc   = "ORC",
    cati  = "CATI",
    hyg   = "Hygiene",
    chlor = "Chlorination",
    lat   = "Latrines",
    vax1  = "Vaccination",
    vax2  = "Vaccination"
  )

  # extract start/end date pairs
  date_rows <- hz_params %>%
    filter(
      grepl("_(start|end)$", parameter),
      !parameter %in% c(
        "outbreak_start", "outbreak_end",
        "first_intervention_date"
      )
    ) %>%
    mutate(
      prefix = sub("_(start|end)$", "", parameter),
      bound  = sub("^.*_", "", parameter)
    ) %>%
    filter(prefix %in% names(type_map), value != "") %>%
    select(hz, prefix, bound, value) %>%
    pivot_wider(names_from = bound, values_from = value) %>%
    mutate(
      start = as.Date(start),
      end   = as.Date(end),
      type  = type_map[prefix]
    ) %>%
    filter(!is.na(start)) # drop interventions with no start date

  # merge vax1 / vax2 rows that share the same type label
  date_rows
}

interventions <- parse_interventions(hz_params)

# Make type a factor with all levels so legends always show every type
all_types <- names(intervention_colours)
interventions$type <- factor(interventions$type, levels = all_types)

# ── load all fit RDS files ───────────────────────────────────────────────
fit_files <- list.files(rds_dir, pattern = "_fit\\.rds$", full.names = TRUE)
# drop the comparative fit
fit_files <- fit_files[!grepl("comparative", fit_files)]

hz_names <- sub("_fit\\.rds$", "", basename(fit_files))

message("Loading ", length(fit_files), " fit files...")

all_fits <- setNames(lapply(fit_files, readRDS), hz_names)

# ── generate fitted trajectories ─────────────────────────────────────────
generate_fit_data <- function(rds, seed = 42, n_draws = 200) {
  fc <- chlaa_forecast_from_fit(
    fit = rds$fit,
    time = rds$observed$time,
    vars = "inc_symptoms_weekly",
    include_cases = TRUE,
    obs_model = "nbinom",
    n_draws = n_draws,
    burnin = 0.25,
    seed = seed,
    dt = 1,
    deterministic = FALSE
  )

  fit_cases <- fc %>%
    filter(variable == "cases") %>%
    left_join(rds$observed %>% select(time, date), by = "time")

  list(fit_cases = fit_cases, observed = rds$observed)
}

message("Generating fitted trajectories...")
all_data <- lapply(all_fits, generate_fit_data)

# ── plotting helpers ─────────────────────────────────────────────────────

# Build intervention segment data for a given HZ
# Nudges start dates that fall within `cluster_window` days of each other
# so dashed lines don't overlap.  Within each cluster the lines are spread
# `spacing` days apart, centred on the cluster midpoint.
get_hz_interventions <- function(hz_name, interventions, y_max,
                                 cluster_window = 10, spacing = 5) {
  intv <- interventions %>% filter(hz == hz_name)
  if (nrow(intv) == 0) {
    return(NULL)
  }
  intv$y_pos <- y_max
  intv <- intv %>% arrange(start, type)

  # greedy clustering: walk sorted starts; any date within

  # `cluster_window` of the previous one joins the same cluster
  n <- nrow(intv)
  cluster_id <- integer(n)
  cl <- 1L
  cluster_id[1] <- cl
  if (n > 1) {
    for (i in 2:n) {
      if (as.integer(intv$start[i] - intv$start[i - 1]) <= cluster_window) {
        cluster_id[i] <- cl
      } else {
        cl <- cl + 1L
        cluster_id[i] <- cl
      }
    }
  }
  intv$cluster <- cluster_id

  # within each cluster, spread lines evenly around the cluster midpoint
  intv <- intv %>%
    group_by(cluster) %>%
    mutate(
      n_cl    = n(),
      idx     = row_number(),
      mid     = mean(start),
      start   = mid + as.difftime((idx - (n_cl + 1) / 2) * spacing,
                                  units = "days")
    ) %>%
    ungroup() %>%
    select(-cluster, -n_cl, -idx, -mid)

  intv
}

# ── Free-text annotation labels for Kirotshe ─────────────────────────────
# Edit this table to add / modify the black arrow+text labels on the
# Kirotshe panel.  Each row draws a horizontal text label with an arrow
# pointing to (arrow_x, arrow_y).
#
#   label    – the text you want displayed
#   x / y    – position of the text anchor (Date / numeric)
#   arrow_x  – x-position the arrow points TO  (typically the intervention
#              start date, i.e. the dashed line)
#   arrow_y  – y-position the arrow points TO
#
# Rows can be added, removed, or reordered freely.

## NOTE: CTC, ORC, and Hygiene all start on the same real-world date in
## Kirotshe (2025-03-24). get_hz_interventions() nudges same-week starts
## +/-5 days apart (in type order: CTC, ORC, Hygiene) so the three dashed
## lines are visible instead of overlapping (landing at 03-19 / 03-24 /
## 03-29) — if `spacing` in get_hz_interventions() ever changes, update
## the arrow_x dates below. Those three labels are stacked to the left
## (left-aligned, hjust = 0) with an arrow to their own line so they don't
## sit on top of one another; CATI/Chlorination keep the original
## centred style (hjust = 0.5).
kirotshe_annotations <- tribble(
  ~label,                   ~x,              ~y,   ~arrow_x,            ~arrow_y,  ~hjust,
  "Cholera Treatment Center (CTC) opened",            "2025-01-08",     195,  "2025-03-19",         180,  0,
  "Oral Rehydration Counter (ORC) opened",             "2025-01-08",     160,  "2025-03-24",         145,  0,
  "Hygiene kits distributed",     "2025-01-08",     125,  "2025-03-29",         110,  0,
  "Localised response (CATI) teams deployed",   "2025-05-20",     200,  "2025-06-23",         180,  0.5,
  "Chlorination points at community water sources",    "2025-04-15",     130,  "2025-05-12",         115,  0.5
) %>%
  mutate(
    x       = as.Date(x),
    arrow_x = as.Date(arrow_x)
  )

# Main (large) plot for Kirotshe with full labels
plot_main <- function(fit_cases, observed, intv, hz_label,
                      annotations = NULL) {
  y_max <- max(c(observed$cases, fit_cases$q0p975), na.rm = TRUE)

  p <- ggplot() +
    # — uncertainty ribbons —
    geom_ribbon(
      data = fit_cases,
      aes(x = date, ymin = q0p025, ymax = q0p975),
      fill = "#6baed6", alpha = 0.25
    ) +
    geom_ribbon(
      data = fit_cases,
      aes(x = date, ymin = q0p25, ymax = q0p75),
      fill = "#6baed6", alpha = 0.45
    ) +
    # — median fit line —
    geom_line(
      data = fit_cases,
      aes(x = date, y = q0p5),
      colour = "#08519c", linewidth = 0.8
    ) +
    # — observed data (grey line) —
    geom_line(
      data = observed,
      aes(x = date, y = cases),
      colour = "grey50", linewidth = 0.4
    )

  # — intervention lines and shaded regions —
  if (!is.null(intv) && nrow(intv) > 0) {
    p <- p +
      geom_vline(
        data = intv,
        aes(xintercept = start, colour = type),
        linetype = "dashed", linewidth = 0.7, alpha = 0.9
      ) +
      geom_rect(
        data = intv %>% filter(!is.na(end)),
        aes(
          xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
          fill = type
        ),
        alpha = 0.04
      )
  }

  # — free-text annotations with arrows (black, horizontal) —
  if (!is.null(annotations) && nrow(annotations) > 0) {
    p <- p +
      geom_segment(
        data = annotations,
        aes(x = x, y = y, xend = arrow_x, yend = arrow_y),
        colour = "black", linewidth = 0.35,
        arrow = arrow(length = unit(0.12, "cm"), type = "closed")
      ) +
      geom_label(
        data = annotations,
        aes(x = x, y = y, label = label, hjust = hjust),
        colour = "black", size = 3.2, fontface = "bold",
        family = "Helvetica",
        fill = "white", linewidth = 0, label.padding = unit(0.15, "lines"),
        vjust = 0.5
      )
  }

  p <- p +
    scale_colour_manual(
      values = intervention_colours,
      limits = names(intervention_colours),
      name   = "Intervention"
    ) +
    scale_fill_manual(
      values = intervention_colours,
      limits = names(intervention_colours),
      name   = "Intervention"
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b %Y",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x     = NULL,
      y     = "Cases per week",
      title = hz_label
    ) +
    theme_minimal(base_size = 13, base_family = "Helvetica") +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.y     = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position  = "none",
      plot.margin      = margin(5, 10, 5, 5)
    )

  p
}

# Small facet-style plot (no text labels, just coloured lines)
plot_small <- function(fit_cases, observed, intv, hz_label) {
  p <- ggplot() +
    geom_ribbon(
      data = fit_cases,
      aes(x = date, ymin = q0p025, ymax = q0p975),
      fill = "#6baed6", alpha = 0.25
    ) +
    geom_ribbon(
      data = fit_cases,
      aes(x = date, ymin = q0p25, ymax = q0p75),
      fill = "#6baed6", alpha = 0.45
    ) +
    geom_line(
      data = fit_cases,
      aes(x = date, y = q0p5),
      colour = "#08519c", linewidth = 0.5
    ) +
    geom_line(
      data = observed,
      aes(x = date, y = cases),
      colour = "grey50", linewidth = 0.3
    )

  if (!is.null(intv) && nrow(intv) > 0) {
    p <- p +
      geom_vline(
        data = intv,
        aes(xintercept = start, colour = type),
        linetype = "dashed", linewidth = 0.55, alpha = 0.9,
        show.legend = FALSE
      ) +
      geom_rect(
        data = intv %>% filter(!is.na(end)),
        aes(
          xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
          fill = type
        ),
        alpha = 0.04,
        show.legend = FALSE
      )
  }

  p <- p +
    scale_colour_manual(values = intervention_colours, drop = FALSE) +
    scale_fill_manual(values = intervention_colours, drop = FALSE) +
    scale_x_date(
      date_breaks = "3 months",
      date_labels = "%b\n%Y",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = NULL, title = hz_label) +
    theme_minimal(base_size = 9, base_family = "Helvetica") +
    theme(
      plot.title        = element_text(face = "bold", size = 9, hjust = 0.5),
      axis.text.x       = element_text(size = 7),
      axis.text.y       = element_text(size = 7),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white", colour = NA),
      legend.position   = "none",
      plot.margin       = margin(2, 4, 2, 4)
    )

  p
}

# ── build all panels ─────────────────────────────────────────────────────

# Title-case helper
hz_title <- function(x) {
  x <- gsub("_", " ", x)
  tools::toTitleCase(x)
}

# Kirotshe main panel
message("Building Kirotshe main panel...")
kiro_data <- all_data[["kirotshe"]]
kiro_intv <- get_hz_interventions(
  "kirotshe", interventions,
  max(
    c(
      kiro_data$observed$cases,
      kiro_data$fit_cases$q0p975
    ),
    na.rm = TRUE
  )
)
p_kiro <- plot_main(
  kiro_data$fit_cases, kiro_data$observed,
  kiro_intv, "Kirotshe",
  annotations = kirotshe_annotations
)

# Other HZs — small panels
other_hz <- setdiff(hz_names, "kirotshe")
other_hz <- sort(other_hz)

message("Building ", length(other_hz), " small panels...")
small_plots <- lapply(other_hz, function(hz) {
  d <- all_data[[hz]]
  intv <- get_hz_interventions(
    hz, interventions,
    max(c(
      d$observed$cases,
      d$fit_cases$q0p975
    ), na.rm = TRUE)
  )
  plot_small(d$fit_cases, d$observed, intv, hz_title(hz))
})

# ── compose with patchwork ───────────────────────────────────────────────
# Layout: Kirotshe on top (spanning full width), small panels below in a grid

n_other <- length(other_hz)
n_cols <- 4
n_rows <- ceiling(n_other / n_cols)

# Combine small plots into a grid
small_grid <- wrap_plots(small_plots, ncol = n_cols) +
  plot_annotation(tag_levels = NULL)

# ── build a visual subtitle legend strip ──────────────────────────────────
# Small ggplot showing line/ribbon symbols next to labels
legend_data <- data.frame(
  xmin  = c(1, 4.5, 9,   14),
  xmax  = c(3, 7.5, 12,  16),
  xmid  = c(2, 6,   10.5, 15),
  label = c("Median", "50% UI", "95% UI", "Observed cases"),
  stringsAsFactors = FALSE
)

p_legend_strip <- ggplot() +
  # 95% UI — medium ribbon + line
  annotate("rect", xmin = 9, xmax = 12, ymin = 0.3, ymax = 0.7,
           fill = "#6baed6", alpha = 0.25) +
  annotate("segment", x = 9, xend = 12, y = 0.5, yend = 0.5,
           colour = "#08519c", linewidth = 0.7) +
  annotate("text", x = 12.2, y = 0.5, label = "95% UI",
           hjust = 0, size = 3.5, family = "Helvetica", colour = "grey30") +
  # 50% UI — darker ribbon + line
  annotate("rect", xmin = 4.5, xmax = 7.5, ymin = 0.3, ymax = 0.7,
           fill = "#6baed6", alpha = 0.45) +
  annotate("segment", x = 4.5, xend = 7.5, y = 0.5, yend = 0.5,
           colour = "#08519c", linewidth = 0.7) +
  annotate("text", x = 7.7, y = 0.5, label = "50% UI",
           hjust = 0, size = 3.5, family = "Helvetica", colour = "grey30") +
  # Median — line only
  annotate("segment", x = 0.5, xend = 3, y = 0.5, yend = 0.5,
           colour = "#08519c", linewidth = 0.8) +
  annotate("text", x = 3.2, y = 0.5, label = "Median",
           hjust = 0, size = 3.5, family = "Helvetica", colour = "grey30") +
  # Observed — grey line
  annotate("segment", x = 15, xend = 17.5, y = 0.5, yend = 0.5,
           colour = "grey50", linewidth = 0.6) +
  annotate("text", x = 17.7, y = 0.5, label = "Observed cases",
           hjust = 0, size = 3.5, family = "Helvetica", colour = "grey30") +
  scale_x_continuous(limits = c(0, 23), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# ── build a manual intervention legend strip ─────────────────────────────
# Each entry: coloured dashed line + shaded rect + label
# Laid out horizontally, centred, matching the intervention_colours palette
intv_legend_entries <- data.frame(
  type  = names(intervention_colours),
  col   = unname(intervention_colours),
  stringsAsFactors = FALSE
)
n_entries <- nrow(intv_legend_entries)
entry_width <- 3.0
gap <- 0.3
total_w <- n_entries * entry_width + (n_entries - 1) * gap
x_start <- (26 - total_w) / 2 + 1.5  # offset right to leave room for title

p_intv_legend <- ggplot()
for (i in seq_len(n_entries)) {
  x0 <- x_start + (i - 1) * (entry_width + gap)
  x1 <- x0 + 1.0  # line/rect width
  x_lab <- x1 + 0.2
  cc <- intv_legend_entries$col[i]
  lab <- intv_legend_entries$type[i]
  # shaded rect
  p_intv_legend <- p_intv_legend +
    annotate("rect", xmin = x0, xmax = x1, ymin = 0.25, ymax = 0.75,
             fill = cc, alpha = 0.25) +
    # dashed line
    annotate("segment", x = x0, xend = x1, y = 0.5, yend = 0.5,
             colour = cc, linewidth = 0.9, linetype = "dashed") +
    # label
    annotate("text", x = x_lab, y = 0.5, label = lab,
             hjust = 0, size = 3.3, family = "Helvetica", colour = "grey20",
             fontface = "bold")
}
p_intv_legend <- p_intv_legend +
  annotate("text", x = x_start - 0.8, y = 0.5, label = "Intervention",
           hjust = 1, size = 3.8, family = "Helvetica", colour = "black",
           fontface = "bold") +
  scale_x_continuous(limits = c(0, 28), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# Stack: subtitle strip, kirotshe, intervention legend, facets
composite <- p_legend_strip / p_kiro / p_intv_legend / small_grid +
  plot_layout(heights = c(0.15, 2, 0.15, n_rows)) +
  plot_annotation(
    title = "Model fits and interventions by health zone",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5,
                                family = "Helvetica")
    )
  )

# ── save ─────────────────────────────────────────────────────────────────
out_path <- file.path(fig_dir, "model_fits_and_interventions.png")

# calculate total height based on number of facet rows
fig_height <- 7 + n_rows * 3 # ~7 in for kirotshe, ~3 in per facet row
fig_width <- 16

message("Saving to: ", out_path)
ggsave(out_path, composite,
  width = fig_width, height = fig_height,
  dpi = 300, bg = "white"
)

# also save as PDF for publication
out_pdf <- file.path(fig_dir, "model_fits_and_interventions.pdf")
message("Saving PDF to: ", out_pdf)
ggsave(out_pdf, composite,
  width = fig_width, height = fig_height,
  device = cairo_pdf, bg = "white"
)

message("Done.")
