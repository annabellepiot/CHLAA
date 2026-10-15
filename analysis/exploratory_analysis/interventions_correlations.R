### add analysis for tuesday. test a github push, will continue to push all for collaboration, but update readme stating "exploratory_analysis" folder is not needed to reproduce our results.

## start the analysis here

# Exploratory analysis of AA intervention indicators vs. IDSR cholera cases (2025).
# One figure per HZ -> figures/groups_<hz>.png, showing cases, case growth and
# each indicator (ordered by intervention group), to see:
#   - whether any intervention is associated with changes in cases / case growth
#   - whether the intervention groups are coherent
#   - whether the model on/off windows (hz_parameters.csv) match reported activity
# One pairwise correlation figure per HZ -> figures/pairs_<hz>.png (Spearman
# heatmap + scatter matrix of cases and indicators).
# Indicator definitions: analysis/data/AA_indicator_codebook.csv
#
# Caveats when reading the figures:
# - orc (treatment counts) follows cases by construction: an association with
#   cases is a sanity check, not evidence of an effect.
# - interventions are often scaled up reactively near the peak; activity that
#   starts after the peak is not evidence of an effect even though cases fall.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(GGally)

data_file <- "Z:/acp25_DIDE/CHLAA/analysis/data/IDSR_AA_indicators.csv"
param_file <- "Z:/acp25_DIDE/CHLAA/analysis/data/hz_parameters.csv"
codebook_file <- "Z:/acp25_DIDE/CHLAA/analysis/data/AA_indicator_codebook.csv"
fig_dir <- "Z:/acp25_DIDE/CHLAA/figures"

# =============================================================================
# Shared setup
# =============================================================================

# ---- weekly data, completed to weeks 1-52 (missing weeks = NA, lines break) ----
week1_date <- as.Date("2024-12-30") # IDSR date of 2025 week 1

dat <- read.csv(data_file, stringsAsFactors = FALSE) %>%
    select(-date) %>%
    complete(hz, week = 1:52) %>%
    mutate(year = 2025, date = week1_date + 7 * (week - 1)) %>%
    arrange(hz, week)

# ---- case growth: log(cases_t + 1) - log(cases_t-1 + 1), centred 3-week mean ----
dat <- dat %>%
    group_by(hz) %>%
    mutate(
        growth_raw = log(cases + 1) - log(lag(cases) + 1),
        growth = as.numeric(stats::filter(growth_raw, rep(1 / 3, 3), sides = 2))
    ) %>%
    ungroup()

# oob_squish: model windows ending in 2026 (Kinshasa) are clipped to the axis
# edge instead of being dropped
x_lim <- range(dat$date) + c(-4, 4)
x_scale <- scale_x_date(
    limits = x_lim, date_breaks = "1 month", date_labels = "%b",
    expand = c(0, 0), oob = scales::oob_squish
)

# ---- intervention groups (match model groups in hz_parameters_definitions.txt) ----
group_ids <- list(
    ctc   = c(2, 4, 13, 14, 24, 25, 26, 37, 42, 43, 45),
    orc   = c(7, 27, 44, 46),
    cati  = c(11, 16, 32, 33),
    hyg   = c(9, 36, 38),
    chlor = c(10, 39),
    lat   = 21
)
group_colours <- c(
    ctc = "#00BFC4", orc = "#FFD700", cati = "#7B68EE",
    hyg = "#FF69B4", chlor = "#32CD32", lat = "#f86d6d",
    other = "#9E9E9E"
)
group_levels <- names(group_colours)

# ---- indicators in long format (id1 = DPS case count, excluded) ----
int_ids <- setdiff(grep("^id\\d+$", names(dat), value = TRUE), "id1")
int_ids <- int_ids[order(as.integer(sub("id", "", int_ids)))]

intervention_groups <- data.frame(id = int_ids, group = "other", stringsAsFactors = FALSE)
for (g in names(group_ids)) intervention_groups$group[intervention_groups$id %in% paste0("id", group_ids[[g]])] <- g
intervention_groups$colour <- group_colours[intervention_groups$group]

int_long <- dat %>%
    select(hz, week, date, all_of(int_ids)) %>%
    pivot_longer(all_of(int_ids), names_to = "id", values_to = "value") %>%
    left_join(intervention_groups, by = "id") %>%
    mutate(
        id = factor(id, levels = int_ids),
        group = factor(group, levels = group_levels)
    )

# keep NA rows (so lines break) but only for ids with any data in that HZ
int_long <- int_long %>%
    group_by(hz, id) %>%
    filter(any(!is.na(value))) %>%
    ungroup()

# ---- top panels: cases (large) + smoothed growth ----
top_panels <- function(h) {
    d <- dat %>% filter(hz == h)

    p_cases <- ggplot(d, aes(date, cases)) +
        geom_col(fill = "grey30", width = 5, na.rm = TRUE) +
        x_scale +
        labs(title = paste0(h, " \u2013 weekly cholera cases (IDSR 2025)"), x = NULL, y = "Cases") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank())

    p_growth <- ggplot(d, aes(date, growth)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        geom_line(colour = "black", na.rm = TRUE) +
        geom_point(colour = "black", size = 0.8, na.rm = TRUE) +
        x_scale +
        labs(x = NULL, y = "Growth\n(3-wk mean)") +
        theme_bw(base_size = 9) +
        theme(panel.grid.minor = element_blank())

    list(cases = p_cases, growth = p_growth)
}

# =============================================================================
# Group figures: cases, growth, indicators by group and model windows -> figures/groups_<hz>.png
# =============================================================================
# Indicators ordered and labelled by group (raw values, group colour). Shaded
# band = model on/off window (<group>_start / <group>_end in hz_parameters.csv).

windows <- read.csv(param_file, stringsAsFactors = FALSE) %>%
    filter(grepl("^(ctc|orc|cati|hyg|chlor|lat)_(start|end)$", parameter)) %>%
    separate(parameter, into = c("group", "bound"), sep = "_") %>%
    mutate(value = as.Date(na_if(value, ""))) %>%
    pivot_wider(names_from = bound, values_from = value) %>%
    filter(!is.na(start), !is.na(end))

for (h in sort(unique(dat$hz))) {
    tp <- top_panels(h)
    int_h <- int_long %>%
        filter(hz == h) %>%
        droplevels()

    if (nrow(int_h) == 0) {
        message("No intervention data for ", h, "; skipping groups figure")
        next
    }

    # one facet per indicator, labelled "group: id" and ordered by group then id
    panel_key <- int_h %>%
        distinct(group, id) %>%
        arrange(group, id) %>%
        mutate(panel = paste0(group, ": ", id))
    panel_key$panel <- factor(panel_key$panel, levels = panel_key$panel)
    int_h <- int_h %>% left_join(panel_key, by = c("group", "id"))

    # model window repeated on every panel of its group
    win_h <- windows %>%
        filter(hz == h) %>%
        inner_join(panel_key %>% mutate(group = as.character(group)),
            by = "group",
            relationship = "many-to-many"
        ) %>%
        mutate(group = factor(group, levels = group_levels))

    n_ids <- nrow(panel_key)

    p_grp <- ggplot(int_h, aes(date, value, colour = group)) +
        geom_rect(
            data = win_h, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = group),
            inherit.aes = FALSE, alpha = 0.15
        ) +
        geom_line(linewidth = 0.5, na.rm = TRUE) +
        geom_point(size = 0.8, na.rm = TRUE) +
        facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "right") +
        scale_colour_manual(values = group_colours, guide = "none") +
        scale_fill_manual(values = group_colours, guide = "none") +
        x_scale +
        labs(
            x = NULL, y = "Indicator value",
            caption = "Shaded band = model intervention window (<group>_start to <group>_end, hz_parameters.csv)"
        ) +
        theme_bw(base_size = 8) +
        theme(
            strip.text.y = element_text(size = 7, angle = 0, margin = margin(1, 2, 1, 2)),
            panel.grid.minor = element_blank(),
            panel.spacing.y = unit(2, "pt")
        )

    p <- tp$cases / tp$growth / p_grp + plot_layout(heights = c(3, 1.2, n_ids * 0.9))

    ggsave(file.path(fig_dir, paste0("groups_", h, ".png")), p,
        width = 12, height = 4.5 + n_ids * 0.9, dpi = 200, limitsize = FALSE
    )
}

# =============================================================================
# Pairwise correlation figures -> figures/pairs_<hz>.png
# =============================================================================
# Adapted from the ggpairs() scatter matrix in the mtcars example, with methods
# suited to short, skewed, autocorrelated weekly series:
# - Spearman rho (rank-based: robust to skew and backlog spikes), pairwise
#   complete weeks, n shown in every tile
# - no p-values / stars: weekly series are autocorrelated (effective n < n) and
#   each figure has many pairs, so p-values would be overconfident
# - linear fit instead of loess (loess overfits with ~8-25 points)
# - case growth is not included: same-week growth added little beyond cases
#   (growth is partly a function of cases_t) and differs from the smoothed
#   growth shown in the groups figures
# - log1p scale for cases and indicators in the scatter panels (display only;
#   Spearman rho is unchanged by monotonic transforms)
# - excluded unit types: rate (unknown reference amounts, mixed scales),
#   no_data, cases (id1 duplicates IDSR cases); * = case_driven_count
# - same-week association only: mostly reflects reactive deployment, not effect

min_weeks <- 8 # minimum weeks an indicator must overlap with cases
min_pairs <- 6 # below this, no rho is shown in a tile

codebook <- read.csv(codebook_file, stringsAsFactors = FALSE, encoding = "UTF-8")
pair_types <- c(
    "weekly_amount", "level", "level_probable", "people_reached",
    "case_driven_count", "rolling_4wk"
)
pair_ids <- intersect(int_ids, codebook$id[codebook$unit_type %in% pair_types])
case_driven <- codebook$id[codebook$unit_type == "case_driven_count"]

upper_spearman <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    ok <- !is.na(x) & !is.na(y)
    n <- sum(ok)
    rho <- if (n >= min_pairs) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA
    lab <- if (is.na(rho)) paste0("n = ", n) else sprintf("ρ = %.2f\nn = %d", rho, n)
    ggplot(data.frame(x = 0, y = 0, rho = rho), aes(x, y)) +
        geom_tile(aes(fill = rho), width = 1, height = 1) +
        annotate("text", x = 0, y = 0, label = lab, size = 2.6) +
        scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            limits = c(-1, 1), na.value = "grey85"
        ) +
        coord_cartesian(expand = FALSE)
}

lower_scatter <- function(data, mapping, ...) {
    ggplot(data, mapping) +
        geom_point(size = 0.8, alpha = 0.7, na.rm = TRUE) +
        geom_smooth(
            method = "lm", formula = y ~ x, se = FALSE,
            linewidth = 0.4, colour = "#B2182B", na.rm = TRUE
        )
}

for (h in sort(unique(dat$hz))) {
    d <- dat %>%
        filter(hz == h, !is.na(cases)) %>%
        select(cases, all_of(pair_ids))

    # indicators with enough weeks and some variation
    keep <- pair_ids[sapply(pair_ids, function(i) {
        v <- d[[i]][!is.na(d[[i]])]
        length(v) >= min_weeks && length(unique(v)) > 1
    })]
    if (length(keep) < 2) {
        message("Fewer than 2 indicators with >= ", min_weeks, " weeks for ", h, "; skipping pairs figure")
        next
    }

    # order by intervention group, then id
    key <- intervention_groups %>%
        filter(id %in% keep) %>%
        mutate(
            group = factor(group, levels = group_levels),
            label = paste0(id, " ", group, ifelse(id %in% case_driven, "*", ""))
        ) %>%
        arrange(group, as.integer(sub("id", "", id)))

    d <- d %>%
        select(cases, all_of(key$id)) %>%
        mutate(across(everything(), log1p))
    n_vars <- ncol(d)

    p <- ggpairs(d,
        columnLabels = c("cases", key$label),
        # na = panels for pairs with no shared weeks (GGally's default na
        # panel errors with ggplot2 >= 4)
        upper = list(continuous = upper_spearman, na = upper_spearman),
        lower = list(continuous = lower_scatter, na = "blank"),
        diag = list(continuous = wrap("barDiag", bins = 8, fill = "grey60", colour = "white")),
        progress = FALSE
    ) +
        labs(
            title = paste0(h, " – pairwise associations (same week, IDSR 2025)"),
            caption = paste(
                "Upper: Spearman ρ (colour, blue = -1, red = +1) and n weeks, pairwise complete;",
                "grey = n <", min_pairs, ". Lower: log1p scale, i.e. log(value + 1), linear fit.",
                "\nNo p-values: weekly series are autocorrelated. Same-week association ≠ effect.",
                "* = case-driven count (correlates with cases by construction)."
            )
        ) +
        theme_bw(base_size = 7) +
        theme(
            strip.text = element_text(size = 6),
            axis.text = element_text(size = 5),
            panel.grid.minor = element_blank()
        )

    size <- max(8, 1.1 * n_vars)
    ggsave(file.path(fig_dir, paste0("pairs_", h, ".png")), p,
        width = size, height = size, dpi = 200
    )
}
